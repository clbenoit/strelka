#!/usr/bin/env python3
#
# Strelka - Small Variant Caller
# Copyright (c) 2009-2018 Illumina, Inc.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#

"""Builder for the fixed Py2/Py3 parity training dataset (change strelka-py3-training).

The canonical training data set (Platinum Genomes NA12878 + ConfidentRegions.bed,
Illumina platgene FTP 2016-1.0) was found unreachable in 2026 (TCP timeout on
192.84.34.14:21, verified). Documented fallback (design.md decision 5): a fixed
synthetic set built from the demo data set + truth, of sufficient size for a
non-degenerate RF. The parity criterion is the equivalence of the two tool
chains, not the clinical quality of the produced model.

Inputs (all from a fixed build of this repo, see README.md):
  - somatic demo run with --reportEVSFeatures  (query VCFs with EVSF features)
  - germline demo run with --reportEVSFeatures, single sample NA12891 (S1)
  - somatic truth VCFs from src/demo/expectedResults
  - demo20.fa reference (5000 bp)

Outputs (scratch/parity/data/):
  - somatic_query_snv.vcf / somatic_query_indel.vcf : real calls + deterministic
    synthetic FP variants (seed 20261010)
  - somatic_truth_snv.vcf / somatic_truth_indel.vcf : demo truth (decompressed)
  - somatic_fp_regions.bed                          : demo20 8 5008 (0-based)
  - germline_query_gvcf.vcf                         : single-sample gVCF with EVSF
  - germline_annotated.vcf                          : synthetic hap.py-style
    two-sample (truth/query) annotated VCF consumed by parseAnnotatedTrainingVcf.py

Determinism: the builder uses only the Python standard library and a fixed RNG
seed; given the same inputs it produces byte-identical outputs. The committed
dataset is pinned by checksums (scratch/parity/checksums.sha256).
"""

import argparse
import gzip
import os
import random

SEED = 20261010
NUC = ["A", "C", "G", "T"]
KMERS = ["AAGC", "TTCG", "GATA", "CTTA", "AGCT", "TTAG", "CAAG", "ATCG"]


def read_ref(fasta, chrom):
    seq = []
    found = False
    with open(fasta) as fid:
        for line in fid:
            if line.startswith(">"):
                found = (chrom in line)
                continue
            if found:
                seq.append(line.strip())
    return "".join(seq)


def read_vcf(path):
    """Return (header_lines, rows) for a (possibly gzipped) VCF."""
    opener = gzip.open if path.endswith(".gz") else open
    header, rows = [], []
    with opener(path, "rt") as fid:
        for line in fid:
            line = line.rstrip("\n")
            if line.startswith("#"):
                header.append(line)
            else:
                rows.append(line.split("\t"))
    return header, rows


def get_info(row, key):
    for field in row[7].split(";"):
        if field.startswith(key + "="):
            return field.split("=", 1)[1]
    return None


def write_vcf(path, header, rows):
    with open(path, "w") as fid:
        for line in header:
            fid.write(line + "\n")
        for row in rows:
            fid.write("\t".join(row) + "\n")


def synth_evf(real_rows, count, rng):
    """Deterministic synthetic EVF rows, data-driven from the observed rows.

    Per feature: if every observed value is integer-valued, draw from the
    distinct observed values with a +-1 jitter (30 %, clamped >= 0); otherwise
    scale an observed value by a uniform factor in [0.4, 1.6] (exact zeros
    stay zero). This preserves the feature space of the real data set.
    """
    nfeat = len(real_rows[0])
    spec = []
    for j in range(nfeat):
        vals = [float(row[j]) for row in real_rows]
        is_int = all(v.is_integer() for v in vals)
        if is_int:
            distinct = sorted({int(v) for v in vals})
            spec.append(("int", distinct))
        else:
            spec.append(("float", vals))
    out = []
    for _ in range(count):
        src = real_rows[rng.randrange(len(real_rows))]
        row = []
        for j in range(nfeat):
            kind, ref = spec[j]
            if kind == "int":
                v = rng.choice(ref)
                if rng.random() < 0.3:
                    v += rng.choice((-1, 1))
                row.append(max(0, v))
            else:
                v = float(src[j])
                if abs(v) < 1e-9:
                    row.append(0)
                else:
                    row.append(round(v * rng.uniform(0.4, 1.6), 5))
        out.append(row)
    return out


def fmt_vals(vals):
    # values may be strings (real VCF fields) or numbers (synthetic rows);
    # normalize to the shortest numeric representation
    out = []
    for v in vals:
        f = float(v)
        out.append(str(int(f)) if f.is_integer() else repr(f))
    return ",".join(out)


def build_somatic(args, rng, outdir):
    ref = read_ref(args.fasta, "demo20")
    assert len(ref) >= 5000, "demo20 too short"

    evsf_dir = args.somatic_evsf_dir
    snv_hdr, snv_rows = read_vcf(os.path.join(evsf_dir, "somatic.snvs.vcf.gz"))
    ind_hdr, ind_rows = read_vcf(os.path.join(evsf_dir, "somatic.indels.vcf.gz"))
    assert all(len(r) == 11 for r in snv_rows + ind_rows)

    # --- SNV query: real calls (8 cols) + synthetic FPs -------------------
    real_snv = [r[:8] for r in snv_rows]
    used = {(r[0], int(r[1]), r[3], r[4]) for r in snv_rows}
    truth_snv = read_vcf(os.path.join(args.expected_dir, "somatic.snvs.vcf.gz"))[1]
    used |= {(r[0], int(r[1]), r[3], r[4]) for r in truth_snv}

    snv_evf_header = [h for h in snv_hdr if "snv_scoring_features=" in h]
    assert len(snv_evf_header) == 1
    nsnv = len(snv_evf_header[0].split("snv_scoring_features=")[1].split(","))
    real_snv_evf = [get_info(r, "EVSF").split(",") for r in snv_rows]
    for evf in real_snv_evf:
        assert len(evf) == nsnv

    snv_synthetic = []
    for i in range(50):
        pos = 100 + 96 * i
        while any((p == pos) for (c, p, a, b) in used):
            pos += 1
        used.add(("demo20", pos, "?", "?"))
        ref_base = ref[pos - 1].upper()
        alt_base = NUC[(NUC.index(ref_base) + 1) % 4]
        qss = rng.randint(3, 40)
        evf = synth_evf(real_snv_evf, 1, rng)[0]
        evf[0] = qss  # SomaticSNVQualityAndHomRefGermlineGenotype tracks QSS
        q = round(rng.uniform(2.0, 12.0), 2)
        dp = rng.randint(5, 40)
        info = ("SOMATIC;QSS=%d;TQSS=1;NT=ref;QSS_NT=%d;TQSS_NT=1;SGT=%s%s->%s%s;"
                "DP=%d;MQ=%s;MQ0=0;ReadPosRankSum=%s;SNVSB=0.00;SomaticEVS=%s;EVSF=%s"
                % (qss, qss, ref_base, ref_base, ref_base, alt_base, dp,
                   evf[3], evf[6], q, fmt_vals(evf)))
        snv_synthetic.append(["demo20", str(pos), ".", ref_base, alt_base,
                              repr(q), "PASS", info])

    snv_rows_out = sorted(
        [list(r) for r in real_snv] + snv_synthetic,
        key=lambda r: (r[0], int(r[1])))
    hdr = (["##fileformat=VCFv4.2",
            "##source=strelka-py3-parity-dataset (builder: tools/build_parity_dataset.py)",
            snv_evf_header[0],
            '##INFO=<ID=EVSF,Number=.,Type=String,Description="EVS scoring features">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"])
    write_vcf(os.path.join(outdir, "somatic_query_snv.vcf"), hdr, snv_rows_out)

    # --- indel query: real calls (8 cols) + synthetic FPs -----------------
    real_ind = [r[:8] for r in ind_rows]
    used_ind = {(r[0], int(r[1]), r[3], r[4]) for r in ind_rows}
    truth_ind = read_vcf(os.path.join(args.expected_dir, "somatic.indels.vcf.gz"))[1]
    used_ind |= {(r[0], int(r[1]), r[3], r[4]) for r in truth_ind}

    ind_evf_header = [h for h in ind_hdr if "indel_scoring_features=" in h]
    assert len(ind_evf_header) == 1
    nind = len(ind_evf_header[0].split("indel_scoring_features=")[1].split(","))
    real_ind_evf = [get_info(r, "EVSF").split(",") for r in ind_rows]
    for evf in real_ind_evf:
        assert len(evf) == nind

    ind_synthetic = []
    for i in range(20):
        pos = 150 + 240 * i
        while any((p == pos) for (c, p, a, b) in used_ind):
            pos += 1
        used_ind.add(("demo20", pos, "?", "?"))
        ref_base = ref[pos - 1].upper()
        kmer = KMERS[i % len(KMERS)]
        if i % 2 == 0:
            ref_a, alt_a = ref_base, ref_base + kmer      # insertion
        else:
            ref_a, alt_a = ref_base + kmer, ref_base      # deletion
        qsi = rng.randint(3, 40)
        evf = synth_evf(real_ind_evf, 1, rng)[0]
        evf[0] = qsi
        q = round(rng.uniform(2.0, 12.0), 2)
        mq = round(float(evf[5]) if abs(float(evf[5])) > 1e-9 else 40.0, 2)
        info = ("SOMATIC;QSI=%d;TQSI=1;NT=ref;QSI_NT=%d;TQSI_NT=1;SGT=ref->het;"
                "MQ=%s;MQ0=0;RU=%s;RC=%d;IC=%d;IHP=%d;SomaticEVS=%s;EVSF=%s"
                % (qsi, qsi, mq, kmer, len(kmer), len(kmer), 1, q, fmt_vals(evf)))
        ind_synthetic.append(["demo20", str(pos), ".", ref_a, alt_a,
                              repr(q), "PASS", info])

    ind_rows_out = sorted(
        [list(r) for r in real_ind] + ind_synthetic,
        key=lambda r: (r[0], int(r[1])))
    hdr = (["##fileformat=VCFv4.2",
            "##source=strelka-py3-parity-dataset (builder: tools/build_parity_dataset.py)",
            ind_evf_header[0],
            '##INFO=<ID=EVSF,Number=.,Type=String,Description="EVS scoring features">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"])
    write_vcf(os.path.join(outdir, "somatic_query_indel.vcf"), hdr, ind_rows_out)

    # --- truth + fp regions -----------------------------------------------
    for src, dst in (("somatic.snvs.vcf.gz", "somatic_truth_snv.vcf"),
                     ("somatic.indels.vcf.gz", "somatic_truth_indel.vcf")):
        hdr, rows = read_vcf(os.path.join(args.expected_dir, src))
        write_vcf(os.path.join(outdir, dst), hdr, rows)
    with open(os.path.join(outdir, "somatic_fp_regions.bed"), "w") as fid:
        fid.write("demo20\t8\t5008\n")

    print("somatic: SNV query %d (real %d + synthetic %d), truth %d; "
          "indel query %d (real %d + synthetic %d), truth %d"
          % (len(snv_rows_out), len(real_snv), 50, len(truth_snv),
             len(ind_rows_out), len(real_ind), 20, len(truth_ind)))


def build_germline(args, rng, outdir):
    ref = read_ref(args.fasta, "demo20")
    evsf_dir = args.germline_evsf_dir
    gvcf_path = os.path.join(evsf_dir, "genome.vcf.gz")
    gvcf_hdr, gvcf_rows = read_vcf(gvcf_path)

    # keep the decompressed gVCF as the filterTrainingVcf input
    write_vcf(os.path.join(outdir, "germline_query_gvcf.vcf"), gvcf_hdr, gvcf_rows)

    snv_feat_line = [h for h in gvcf_hdr if "snv_scoring_features=" in h][0]
    ind_feat_line = [h for h in gvcf_hdr if "indel_scoring_features=" in h][0]
    nsnv = len(snv_feat_line.split("snv_scoring_features=")[1].split(","))
    nind = len(ind_feat_line.split("indel_scoring_features=")[1].split(","))

    # real called variants (single sample S1)
    calls = []
    for row in gvcf_rows:
        if not any(g in row[9] for g in ("0/1", "1/0", "0/2")):
            continue
        evsf = get_info(row, "EVSF")
        assert evsf is not None
        is_indel = ("CIGAR=" in row[7])
        vals = evsf.split(",")
        assert len(vals) == (nind if is_indel else nsnv), (row[1], len(vals))
        calls.append({"chrom": row[0], "pos": int(row[1]), "ref": row[3],
                      "alt": row[4], "qual": row[5], "evsf": vals,
                      "type": "INDEL" if is_indel else "SNP"})
    assert calls, "no germline calls found"
    used = {(c["chrom"], c["pos"], c["ref"], c["alt"]) for c in calls}
    real_evf = {"SNP": [c["evsf"] for c in calls if c["type"] == "SNP"],
                "INDEL": [c["evsf"] for c in calls if c["type"] == "INDEL"]}
    assert real_evf["SNP"] and real_evf["INDEL"]

    # synthetic FPs: 40 SNV + 10 INDEL at deterministic positions
    synthetic = []
    for i in range(50):
        pos = 200 + 90 * i
        while any((p == pos) for (c, p, a, b) in used):
            pos += 1
        used.add(("demo20", pos, "?", "?"))
        is_indel = (i % 5 == 4)
        ref_base = ref[pos - 1].upper()
        if is_indel:
            kmer = KMERS[i % len(KMERS)]
            ref_a, alt_a = (ref_base + kmer, ref_base) if i % 2 else (ref_base, ref_base + kmer)
            vtype = "INDEL"
        else:
            alt_base = NUC[(NUC.index(ref_base) + 1) % 4]
            ref_a, alt_a = ref_base, alt_base
            vtype = "SNP"
        evf = synth_evf(real_evf[vtype], 1, rng)[0]
        q = round(rng.uniform(2.0, 12.0), 2)
        synthetic.append({"chrom": "demo20", "pos": pos, "ref": ref_a,
                          "alt": alt_a, "qual": repr(q), "evsf": evf,
                          "type": vtype})

    # FN: shifted copies of every 5th real call (truth-only, no EVSF)
    fn_rows = []
    for idx, c in enumerate(calls):
        if idx % 5 == 0:
            assert all(p != c["pos"] + 1 for (ch, p, a, b) in used)
            fn_rows.append({"chrom": c["chrom"], "pos": c["pos"] + 1,
                            "ref": c["ref"], "alt": c["alt"],
                            "type": c["type"]})
            used.add((c["chrom"], c["pos"] + 1, c["ref"], c["alt"]))

    # two-sample hap.py-style annotated VCF: FORMAT BD:BVT, samples truth/query
    def cigar(ref_a, alt_a):
        diff = abs(len(alt_a) - len(ref_a))
        if diff == 0:
            return ""
        return "1M%d%s" % (diff, "I" if len(alt_a) > len(ref_a) else "D")

    lines = []
    for c in calls:  # TP rows
        info = "EVSF=" + fmt_vals(c["evsf"])
        cg = cigar(c["ref"], c["alt"])
        if cg:
            info += ";CIGAR=" + cg
        lines.append(["demo20", str(c["pos"]), ".", c["ref"], c["alt"],
                      c["qual"], "PASS", info, "BD:BVT",
                      "TP:" + c["type"], "TP:" + c["type"]])
    for c in synthetic:  # FP rows
        info = "EVSF=" + fmt_vals(c["evsf"])
        cg = cigar(c["ref"], c["alt"])
        if cg:
            info += ";CIGAR=" + cg
        lines.append(["demo20", str(c["pos"]), ".", c["ref"], c["alt"],
                      c["qual"], "PASS", info, "BD:BVT",
                      "FP:" + c["type"], "FP:" + c["type"]])
    for c in fn_rows:  # FN rows: truth only, no EVSF
        lines.append(["demo20", str(c["pos"]), ".", c["ref"], c["alt"],
                      ".", "PASS", ".", "BD:BVT",
                      "FN:" + c["type"], "NOCALL:" + c["type"]])
    lines.sort(key=lambda r: int(r[1]))

    hdr = ["##fileformat=VCFv4.2",
           "##source=strelka-py3-parity-dataset (synthetic hap.py-style annotation)",
           "##FORMAT=<ID=BD,Number=1,Type=String,Description=\"hap.py-style label (TP/FP/FN)\">",
           "##FORMAT=<ID=BVT,Number=1,Type=String,Description=\"variant type (SNP/INDEL)\">",
           "##INFO=<ID=EVSF,Number=.,Type=String,Description=\"EVS scoring features\">",
           snv_feat_line,
           ind_feat_line,
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttruth\tquery"]
    write_vcf(os.path.join(outdir, "germline_annotated.vcf"), hdr, lines)

    print("germline: gVCF %d rows, calls %d (TP), synthetic FP %d, FN %d; "
          "annotated VCF %d rows"
          % (len(gvcf_rows), len(calls), len(synthetic), len(fn_rows), len(lines)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--somatic-evsf-dir", required=True)
    parser.add_argument("--germline-evsf-dir", required=True)
    parser.add_argument("--expected-dir", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    rng = random.Random(SEED)
    os.makedirs(args.out, exist_ok=True)
    build_somatic(args, rng, args.out)
    build_germline(args, rng, args.out)
    print("dataset written to", args.out)


if __name__ == "__main__":
    main()

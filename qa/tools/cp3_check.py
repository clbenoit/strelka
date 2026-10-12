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
# but WITHOUT OTHER WARRANTY; without even the implied warranty of
# MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#


"""Checkpoint 3 semantic check (change strelka-qa-noregression).

For each VCF pair (py2-model run vs py3-model run), verify the non-regression
criterion used at parity validation (scratch/parity/report.md, CP3):

  - variant SET parity: identical (CHROM, POS, REF, ALT) — the pass criterion
  - score deltas: per variant, QUAL and every numeric FORMAT field
    (SomaticEVS, GQ, GQX, ...) compared as floats — quantified, informational

Exit code 0 iff every pair has identical variant sets.

Usage:
    cp3_check.py --somatic-snvs A B --somatic-indels A B --germline A B
"""

import argparse
import gzip
import sys


def read_variants(path):
    """Return { (chrom,pos,ref,alt): (qual, [format_values_per_sample]) }."""
    variants = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fid:
        header = []
        for line in fid:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line.strip().split("\t")
                continue
            fields = line.rstrip("\n").split("\t")
            key = (fields[0], fields[1], fields[3], fields[4])
            qual = fields[5]
            fmt_idx = header.index("FORMAT") if "FORMAT" in header else 8
            samples = fields[fmt_idx + 1:]
            variants[key] = (qual, samples)
    return variants


def numeric_values(sample_values):
    """Extract the numeric (non-GT) fields of per-sample FORMAT values."""
    out = []
    for sample in sample_values:
        parts = sample.split(":")
        # parts[0] is GT (or the first non-numeric field); compare the rest
        for part in parts[1:]:
            part = part.split(",")[0]
            try:
                out.append(float(part))
            except ValueError:
                pass
    return out


def compare_pair(label, path_a, path_b):
    va = read_variants(path_a)
    vb = read_variants(path_b)
    only_a = set(va) - set(vb)
    only_b = set(vb) - set(va)
    common = set(va) & set(vb)
    if only_a or only_b:
        detail = []
        if only_a:
            detail.append("only A: %s" % sorted(only_a)[:5])
        if only_b:
            detail.append("only B: %s" % sorted(only_b)[:5])
        return ("SET-DIFF (%d/%d vs %d/%d; %s)"
                % (len(only_a), len(va), len(only_b), len(vb),
                   "; ".join(detail)), False)

    changed = 0
    max_delta = 0.0
    for key in common:
        qual_a, samp_a = va[key]
        qual_b, samp_b = vb[key]
        try:
            delta = abs(float(qual_a) - float(qual_b))
            if delta > 0:
                changed += 1
            max_delta = max(max_delta, delta)
        except ValueError:
            pass
        nums_a = numeric_values(samp_a)
        nums_b = numeric_values(samp_b)
        for x, y in zip(nums_a, nums_b):
            d = abs(x - y)
            if d > 0:
                changed += 1
            max_delta = max(max_delta, d)
    return ("%d/%d variants identical; %d score fields changed "
            "(max |d| %.2f)" % (len(common), len(va), changed, max_delta),
            True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--somatic-snvs", nargs=2, required=True)
    parser.add_argument("--somatic-indels", nargs=2, required=True)
    parser.add_argument("--germline", nargs=2, required=True)
    args = parser.parse_args()

    pairs = [
        ("somatic.snvs", args.somatic_snvs),
        ("somatic.indels", args.somatic_indels),
        ("germline", args.germline),
    ]
    ok = True
    parts = []
    for label, (a, b) in pairs:
        result, pair_ok = compare_pair(label, a, b)
        ok = ok and pair_ok
        parts.append("%s: %s" % (label, result))
    print(" ; ".join(parts))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()

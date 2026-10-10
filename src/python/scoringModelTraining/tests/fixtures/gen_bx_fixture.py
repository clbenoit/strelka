# -*- coding: utf-8 -*-
"""Throwaway fixture generator (run under the py2 reference toolchain).

Builds a BedIntervalTree with the UNCONVERTED (py2/HEAD) code from the main
checkout, using the exact consumer pattern of vcf_to_feature_csv.py:

    labeledIntervals.addFromBed(fp_bed, "FP")
    labeledIntervals.addFromBed(ambig_bed, lambda xe: xe[4])

then answers every query from bedintervaltree_parity_queries.txt and writes
the expected results (sorted, order-independent) to
bedintervaltree_parity_bx_results.json. The py3 unit test
(test_intervaltree_parity.py) compares the in-repo implementation against
this fixture; no py2 dependency is needed by the py3 suite.
"""
import json
import os
import sys

FIX = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, "/workspace/repos/strelka/src/python/scoringModelTraining/somatic/lib")

import evs
from evs.tools.bedintervaltree import BedIntervalTree  # noqa: E402

tree = BedIntervalTree()
tree.addFromBed(os.path.join(FIX, "bedintervaltree_parity_fp.bed"), "FP")
tree.addFromBed(os.path.join(FIX, "bedintervaltree_parity_ambig.bed"),
                lambda xe: xe[4])

counts = {
    "total": tree.count(),
    "FP": tree.count("FP"),
    "repeat": tree.count("repeat"),
    "lowqual": tree.count("lowqual"),
}

queries = []
with open(os.path.join(FIX, "bedintervaltree_parity_queries.txt")) as fid:
    for line in fid:
        line = line.strip()
        if not line:
            continue
        chrom, start, end = line.split("\t")
        start, end = int(start), int(end)
        res = tree.intersect(chrom, start, end)
        rows = sorted(
            [i.start, i.end, [unicode(v) for v in i.value], unicode(i.chrom)]
            for i in res
        )
        queries.append({"chrom": chrom, "start": start, "end": end, "results": rows})

out = os.path.join(FIX, "bedintervaltree_parity_bx_results.json")
with open(out, "w") as fid:
    json.dump({"counts": counts, "queries": queries}, fid, indent=1)
    fid.write("\n")

print("wrote %s" % out)
print("counts:", counts)
for q in queries:
    print("%s %d %d -> %d" % (q["chrom"], q["start"], q["end"], len(q["results"])))

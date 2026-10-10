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
#

"""Checkpoint 2: compare py2 vs py3 exported EVS model JSONs.

Checks, per model:
  - metadata parity (Features list, ModelType, FilterCutoff, Calibration)
  - number of trees
  - tree structure parity (children layout, split features) and quantified
    divergence (node counts, threshold agreement where structure matches)
  - node_votes semantics: py2 (sklearn 0.19.2) = class counts per node,
    py3 (sklearn 1.9.1) = normalized class probabilities (sum 1.0). The C++
    RandomForestModel consumes the ratio vote.left/(vote.left+vote.right),
    which is scale-invariant, so both layouts are drop-in equivalent.

Usage: compare_models.py <model-spec> ...
  <model-spec> = <tag2>:<tag3>:<path2>:<path3>
    e.g. py2:py3:results/py2/model_somatic_snv.json:results/py3/model_somatic_snv.json
"""

import json
import sys


def load_model(path):
    with open(path) as fid:
        data = json.load(fid)
    domain = list(data["CalibrationModels"].keys())[0]
    section = list(data["CalibrationModels"][domain].keys())[0]
    return data["CalibrationModels"][domain][section], domain, section


def tree_structure(model):
    out = []
    for t in model["Model"]:
        nodes = len(t["tree"])
        leaves = sum(1 for v in t["tree"].values() if v[0] == -1)
        feats = sorted({v[0] for v in t["decisions"].values() if v[0] != -1})
        out.append({"nodes": nodes, "leaves": leaves, "features": feats,
                    "tree": t["tree"], "decisions": t["decisions"]})
    return out


def compare(name, path2, path3):
    m2, domain2, section2 = load_model(path2)
    m3, domain3, section3 = load_model(path3)
    if (domain2, section2) != (domain3, section3):
        print("  WARNING: model sections differ: %s/%s vs %s/%s"
              % (domain2, section2, domain3, section3))

    print("== %s (%s/%s)" % (name, domain2, section2))
    ok = True
    for key in ("Features", "ModelType", "FilterCutoff", "Calibration"):
        same = (m2[key] == m3[key])
        ok &= same
        print("  meta %-13s : %s" % (key, "IDENTICAL" if same else "DIFFERS: %s vs %s" % (m2[key], m3[key])))
    n2, n3 = len(m2["Model"]), len(m3["Model"])
    ok &= (n2 == n3)
    print("  n_trees         : %d vs %d %s" % (n2, n3, "" if n2 == n3 else "(DIFFERS)"))

    t2, t3 = tree_structure(m2), tree_structure(m3)
    same_struct = same_thresh = 0
    node_counts = []
    for a, b in zip(t2, t3):
        node_counts.append((a["nodes"], b["nodes"]))
        if a["tree"] == b["tree"] and all(
                (x[0] == y[0]) for x, y in zip(
                    (v for v in a["decisions"].values()), (v for v in b["decisions"].values()))):
            same_struct += 1
            if all(abs(x[1] - y[1]) < 1e-6 for x, y in zip(
                    (v for v in a["decisions"].values()), (v for v in b["decisions"].values()))):
                same_thresh += 1
    print("  trees structurally identical (same children + split features) : %d/%d"
          % (same_struct, n2))
    print("  trees identical incl. thresholds (< 1e-6)                     : %d/%d"
          % (same_thresh, n2))
    if node_counts:
        deltas = [abs(a - b) for a, b in node_counts]
        print("  per-tree node count: py2 %d..%d, py3 %d..%d, max |delta| %d"
              % (min(a for a, _ in node_counts), max(a for a, _ in node_counts),
                 min(b for _, b in node_counts), max(b for _, b in node_counts),
                 max(deltas)))
        same_nodes = sum(1 for d in deltas if d == 0)
        print("  trees with identical node count                             : %d/%d" % (same_nodes, n2))

    # node_votes semantics
    def vote_stats(model):
        vals = [v for t in model["Model"] for v in t["node_votes"].values()]
        sums = [sum(v) for v in vals]
        ints = all(float(v[0]).is_integer() and float(v[1]).is_integer() for v in vals)
        return len(vals), min(sums), max(sums), ints

    n2v, lo2, hi2, ints2 = vote_stats(m2)
    n3v, lo3, hi3, ints3 = vote_stats(m3)
    print("  node_votes py2: %d entries, sum range [%.4f, %.4f], all-integer=%s (counts expected)"
          % (n2v, lo2, hi2, ints2))
    print("  node_votes py3: %d entries, sum range [%.4f, %.4f], all-integer=%s (probabilities expected)"
          % (n3v, lo3, hi3, ints3))
    root2 = m2["Model"][0]["node_votes"]["0"]
    root3 = m3["Model"][0]["node_votes"]["0"]
    print("  tree0 root votes: py2 %s | py3 %s (scale-invariant ratio for C++: %.6f vs %.6f)"
          % (root2, root3,
             root2[0] / sum(root2), root3[0] / sum(root3)))
    return ok


def main():
    for spec in sys.argv[1:]:
        tag2, tag3, path2, path3 = spec.split(":")
        compare("%s vs %s" % (tag2, tag3), path2, path3)
        print()


if __name__ == "__main__":
    main()

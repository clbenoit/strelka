#!/usr/bin/env python
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

"""Checkpoint 2b: score the same labeled CSV with each tool chain's own
exported pickle and emit one line per variant (POS, REF, ALT, ptag, P(TP)).

The pickles are not interchangeable across sklearn versions, so each
tool chain scores the CSV with the model it trained itself. The emitted
lines are compared across tool chains: ptag agreement and the P(TP) delta
quantify the functional impact of the tree-level divergences reported by
compare_models.py.

Usage:
    score_check.py --libdir <lib> --pkl FILE --csv FILE --features NAME
"""

import argparse
import os
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--libdir", required=True)
    parser.add_argument("--pkl", required=True)
    parser.add_argument("--csv", required=True)
    parser.add_argument("--features", required=True,
                        help="FeatureSet name (somatic.snv / germline.snv / germline.indel)")
    args = parser.parse_args()

    sys.path.insert(0, os.path.abspath(args.libdir))
    import evs  # noqa: E402
    import evs.features  # noqa: E402
    import pandas  # noqa: E402

    fset = evs.features.FeatureSet.make(args.features)
    features = fset.trainingfeatures()

    model = evs.EVSModel.createFromFile(args.pkl)
    df = pandas.read_csv(args.csv)
    df = df[df["tag"] != "FN"]
    has_alleles = ("REF" in df.columns) and ("ALT" in df.columns)
    res = model.classify(df, features)
    for _, row in res.iterrows():
        if has_alleles:
            print("%s\t%s\t%s\t%s\t%.6f"
                  % (row["POS"], row["REF"], row["ALT"], row["ptag"], row["qual"]))
        else:
            print("%s\t-\t-\t%s\t%.6f"
                  % (row["POS"], row["ptag"], row["qual"]))


if __name__ == "__main__":
    main()

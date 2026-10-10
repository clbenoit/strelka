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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#

"""Deterministic post-processing of the somatic indel feature CSV.

``vcf_to_feature_csv.py`` emits CHROM, POS, REF, ALT, <EVSF features>, tag,
but ``StrelkaRFIndel.train`` additionally consumes two derived columns:

  - NT        : "ref" for every row (the demo data set is all NT=ref)
  - INDELTYPE : 1 = insertion (len(ALT) > len(REF)), 2 = deletion
                (len(ALT) < len(REF))

Both columns are derived purely from REF/ALT; given the same input CSV this
script produces a byte-identical output, so it is run once and the result is
shared by the py2 and py3 reference runs (parity is checked on the upstream
CSV emitted by vcf_to_feature_csv.py).
"""

import csv
import sys


def main():
    if len(sys.argv) != 3:
        sys.exit("usage: postprocess_indel_csv.py <in.csv> <out.csv>")
    with open(sys.argv[1]) as fin:
        reader = csv.reader(fin)
        header = next(reader)
        if "NT" in header or "INDELTYPE" in header:
            sys.exit("input already has NT/INDELTYPE columns")
        try:
            iref, ialt = header.index("REF"), header.index("ALT")
        except ValueError:
            sys.exit("REF/ALT columns not found")
        rows = list(reader)
    out = [header + ["NT", "INDELTYPE"]]
    for row in rows:
        itype = "1" if len(row[ialt]) > len(row[iref]) else "2"
        out.append(row + ["ref", itype])
    with open(sys.argv[2], "w") as fout:
        csv.writer(fout, lineterminator="\n").writerows(out)
    print("wrote %s (%d rows)" % (sys.argv[2], len(rows)))


if __name__ == "__main__":
    main()

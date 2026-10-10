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

"""bx-intervals parity test for the in-repo Python 3 interval tree.

The fixture (``fixtures/bedintervaltree_parity_bx_results.json``) was generated
by the Python 2 reference toolchain (bx-python 0.8.7) from the unconverted
code, using the exact consumer pattern of ``vcf_to_feature_csv.py``::

    labeledIntervals.addFromBed(fp_bed, "FP")
    labeledIntervals.addFromBed(ambig_bed, lambda xe: xe[4])

The test runs the same inputs through the in-repo replacement
(``evs.tools.intervaltree`` + ``evs.tools.bedintervaltree``) under Python 3
and asserts that every query returns the same set of intervals (start, end,
value, chrom). The result *order* of ``find`` is not a consumed contract
(callers aggregate values into sets/flags), so comparison is on sorted rows.

No Python 2 dependency is needed to run this suite.
"""

import json
import os
import sys
import unittest

_TESTSDIR = os.path.dirname(os.path.abspath(__file__))
_FIXDIR = os.path.join(_TESTSDIR, "fixtures")
# the somatic tree carries the consumer (vcf_to_feature_csv.py); both trees
# (somatic/germline) ship an identical copy of the implementation
sys.path.insert(0, os.path.join(_TESTSDIR, "..", "somatic", "lib"))


class TestBedIntervalTreeBxParity(unittest.TestCase):
    """In-repo replacement vs bx-python 0.8.7 reference fixture."""

    @classmethod
    def setUpClass(cls):
        from evs.tools.bedintervaltree import BedIntervalTree

        cls.BedIntervalTree = BedIntervalTree
        with open(os.path.join(_FIXDIR, "bedintervaltree_parity_bx_results.json")) as fid:
            cls.fixture = json.load(fid)
        cls.queries = []
        with open(os.path.join(_FIXDIR, "bedintervaltree_parity_queries.txt")) as fid:
            for line in fid:
                line = line.strip()
                if line:
                    chrom, start, end = line.split("\t")
                    cls.queries.append((chrom, int(start), int(end)))

    def _buildTree(self):
        tree = self.BedIntervalTree()
        tree.addFromBed(os.path.join(_FIXDIR, "bedintervaltree_parity_fp.bed"), "FP")
        tree.addFromBed(os.path.join(_FIXDIR, "bedintervaltree_parity_ambig.bed"),
                        lambda xe: xe[4])
        return tree

    def test_query_set_parity(self):
        """Every query returns the same set of (start, end, value, chrom)."""
        tree = self._buildTree()
        by_query = {(q["chrom"], q["start"], q["end"]): q["results"]
                    for q in self.fixture["queries"]}
        self.assertEqual(len(by_query), len(self.queries))
        for chrom, start, end in self.queries:
            def canon(rows):
                return sorted((int(s), int(e), tuple(v), str(c)) for s, e, v, c in rows)

            expected = canon(by_query[(chrom, start, end)])
            got = tree.intersect(chrom, start, end)
            rows = canon((i.start, i.end, list(i.value), i.chrom) for i in got)
            self.assertEqual(
                rows, expected,
                "query %s %d %d: got %s, expected %s" % (chrom, start, end, rows, expected))

    def test_counts_parity(self):
        """Interval counts per label match the reference (duplicates included)."""
        tree = self._buildTree()
        counts = self.fixture["counts"]
        self.assertEqual(tree.count(), counts["total"])
        for label in ("FP", "repeat", "lowqual"):
            self.assertEqual(tree.count(label), counts[label])

    def test_value_and_chrom_preserved(self):
        """value = [label] + bed columns 3+, chrom attribute kept (consumer contract)."""
        tree = self._buildTree()
        # FP entries: 3-column BED, fixed label -> value == ["FP"]
        # (chr1 5..10: only the FP interval [1, 51) overlaps)
        res = tree.intersect("chr1", 5, 10)
        self.assertEqual(len(res), 1)
        self.assertEqual(res[0].value, ["FP"])
        self.assertEqual(res[0].chrom, "chr1")
        # ambiguous entries: 5-column BED, label = column 4 -> value = [label] + cols 3..
        # (chr1 55..60: only the ambiguous interval [41, 61) overlaps)
        res = tree.intersect("chr1", 55, 60)
        self.assertEqual(len(res), 1)
        self.assertEqual(res[0].value, ["lowqual", "editing", "lowqual"])
        self.assertEqual(res[0].chrom, "chr1")


class TestIntervalTreeBoundaries(unittest.TestCase):
    """Half-open [start, end) semantics of the raw in-repo IntervalTree.

    Mirrors the boundary probes run against bx-python 0.8.7 (see change
    strelka-py3-training): touching boundaries do not match, strict point
    containment for empty queries.
    """

    @classmethod
    def setUpClass(cls):
        from evs.tools.intervaltree import Interval, IntervalTree

        cls.Interval = Interval
        cls.IntervalTree = IntervalTree

    def _tree(self):
        t = self.IntervalTree()
        for start, end, value in ((10, 20, "a"), (30, 40, "b"), (50, 60, "c")):
            t.add_interval(self.Interval(start, end, value=value, chrom="x"))
        return t

    def _found(self, t, start, end):
        return sorted((i.start, i.end, i.value) for i in t.find(start, end))

    def test_touching_end_does_not_match(self):
        self.assertEqual(self._found(self._tree(), 20, 25), [])

    def test_touching_start_does_not_match(self):
        self.assertEqual(self._found(self._tree(), 5, 10), [])

    def test_overlapping_matches(self):
        self.assertEqual(self._found(self._tree(), 15, 25), [(10, 20, "a")])

    def test_empty_query_point_containment(self):
        t = self.IntervalTree()
        t.add_interval(self.Interval(5, 15, value="a", chrom="c"))
        t.add_interval(self.Interval(15, 25, value="b", chrom="c"))
        t.add_interval(self.Interval(0, 100, value="c", chrom="c"))
        # point 15: contained only by [0, 100) (strict containment)
        self.assertEqual(self._found(t, 15, 15), [(0, 100, "c")])

    def test_duplicates_all_returned(self):
        t = self.IntervalTree()
        t.add_interval(self.Interval(10, 20, value="a", chrom="c"))
        t.add_interval(self.Interval(10, 20, value="b", chrom="c"))
        self.assertEqual(self._found(t, 5, 30),
                         [(10, 20, "a"), (10, 20, "b")])

    def test_result_count_large_bed(self):
        """sanity: 10k intervals, one wide query returns all overlapping."""
        import random
        t = self.IntervalTree()
        rng = random.Random(1234)
        n = 10000
        for i in range(n):
            start = rng.randrange(0, 100000)
            t.add_interval(self.Interval(start, start + rng.randrange(1, 50),
                                         value=i, chrom="c"))
        # [0, 100000): every interval with start < 100000 and end > 0
        got = t.find(0, 100000)
        self.assertGreater(len(got), n * 0.9)


class TestSharedToolFilesIdentical(unittest.TestCase):
    """The somatic and germline trees ship byte-identical copies of the tool
    modules (as in HEAD); guard against drift (the germline tree's runtime
    behavior is therefore covered by the somatic tests above)."""

    def test_tool_files_byte_identical(self):
        for name in ("intervaltree.py", "bedintervaltree.py", "io.py", "vcf.py"):
            somatic = os.path.join(_TESTSDIR, "..", "somatic", "lib", "evs", "tools", name)
            germline = os.path.join(_TESTSDIR, "..", "germline", "lib", "evs", "tools", name)
            with open(somatic, "rb") as fs, open(germline, "rb") as fg:
                self.assertEqual(fs.read(), fg.read(), name)


if __name__ == "__main__":
    unittest.main()

#
# Strelka - Small Variant Caller
# Copyright (c) 2009-2018 Illumina, Inc.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option) any later version.
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

"""
Minimal Python 3 in-repo replacement for the small surface of
``bx.intervals.intersection`` consumed by ``bedintervaltree.BedIntervalTree``:

- ``Interval(start, end, value=None, chrom=None)`` with attributes
  ``.start``, ``.end``, ``.value``, ``.chrom``;
- ``IntervalTree.add_interval(interval)``;
- ``IntervalTree.find(start, end)``.

Overlap semantics are identical to bx-python 0.8.7 (half-open intervals):
``find(start, end)`` returns every interval ``i`` with ``i.start < end`` and
``i.end > start``.

The result order of ``find`` is a traversal order and is *not* part of the
consumed contract (the single caller aggregates interval values into sets
and flags); parity with the bx-python reference is asserted as a sorted
multiset by the unit tests under ``somatic/lib/evs/tools/tests``.

The structure is a randomized binary search tree (treap) keyed by
``(start, end, insertion sequence)`` and augmented with the maximum ``end``
in each subtree, giving O(log n + k) per query in expectation. The random
priorities are drawn from a fixed-seed generator so identical insertion
sequences produce identical trees (reproducible runs).
"""

import random


class Interval(object):
    """An interval [start, end) carrying an arbitrary value and an optional chromosome.

    Mirrors the bx.intervals.intersection.Interval constructor surface
    (``Interval(start, end, value=..., chrom=...)``).
    """

    __slots__ = ("start", "end", "value", "chrom")

    def __init__(self, start, end, value=None, chrom=None):
        self.start = start
        self.end = end
        self.value = value
        self.chrom = chrom

    def __eq__(self, other):
        return (isinstance(other, Interval)
                and self.start == other.start
                and self.end == other.end
                and self.value == other.value
                and self.chrom == other.chrom)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash((self.start, self.end))

    def __repr__(self):
        return "Interval(start=%r, end=%r, value=%r, chrom=%r)" % (
            self.start, self.end, self.value, self.chrom)


class _Node(object):
    __slots__ = ("interval", "key", "priority", "left", "right", "max_end")


def _max_end(node):
    return node.max_end if node is not None else None


class IntervalTree(object):
    """Augmented treap over Interval, supporting overlap queries.

    Parity with bx-python 0.8.7:
    - ``add_interval`` accepts any Interval (duplicates allowed);
    - ``find(start, end)`` returns all intervals overlapping the half-open
      query range [start, end), i.e. ``i.start < end and i.end > start``.
    """

    # fixed seed: same insertion sequence -> same tree -> reproducible results
    _RNG_SEED = 0

    def __init__(self):
        self._root = None
        self._seq = 0
        self._rng = random.Random(self._RNG_SEED)

    def add_interval(self, interval):
        """Add an interval to the tree (bx-compatible method name)."""
        node = _Node()
        node.interval = interval
        node.key = (interval.start, interval.end, self._seq)
        self._seq += 1
        node.priority = self._rng.random()
        node.left = None
        node.right = None
        node.max_end = interval.end

        self._root = self._insert(self._root, node)
        return interval

    # -- treap machinery ----------------------------------------------------

    def _insert(self, node, new):
        if node is None:
            return new
        if new.key < node.key:
            node.left = self._insert(node.left, new)
            if node.left.priority > node.priority:
                node = self._rotate_right(node)
        else:
            node.right = self._insert(node.right, new)
            if node.right.priority > node.priority:
                node = self._rotate_left(node)
        node.max_end = node.interval.end
        left_end = _max_end(node.left)
        if left_end is not None and left_end > node.max_end:
            node.max_end = left_end
        right_end = _max_end(node.right)
        if right_end is not None and right_end > node.max_end:
            node.max_end = right_end
        return node

    @staticmethod
    def _rotate_right(p):
        q = p.left
        p.left = q.right
        q.right = p
        return q

    @staticmethod
    def _rotate_left(p):
        q = p.right
        p.right = q.left
        q.left = p
        return q

    # -- query --------------------------------------------------------------

    def find(self, start, end):
        """Return all intervals overlapping the half-open range [start, end).

        Overlap rule (identical to bx-python 0.8.7):
        ``interval.start < end and interval.end > start``.
        """
        out = []
        self._collect(self._root, start, end, out)
        return out

    def _collect(self, node, start, end, out):
        # prune: no interval in this subtree ends after the query starts
        if node is None or node.max_end <= start:
            return
        self._collect(node.left, start, end, out)
        if node.interval.start < end and node.interval.end > start:
            out.append(node.interval)
        # every interval in the right subtree starts at or after node.interval.start
        if node.interval.start < end:
            self._collect(node.right, start, end, out)

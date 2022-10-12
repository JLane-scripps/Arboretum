from bisect import bisect, bisect_left
from collections import deque
from dataclasses import dataclass, field
from typing import List

from boundary import Boundary, psm_attributes_in_bound
from forest.psmtree import PsmTree
from psm import PSM


@dataclass
class PsmSortedList(PsmTree):
    """
    PsmSortedList can use either binary search bisect functions or a deque to be much faster than PsmList...
    but still long. Bisect had significantly shorter search times but horrendous add times
    Deque had horrendous search times (after linked lists of 100k psm's) but slightly faster add times.
    """
    tree: List[PSM] = field(default_factory=lambda: list())  # The tree is composed of lists of PSM's
    mz_list: [float] = field(default_factory=lambda: list())  # this is a list of just mz values, for matching indexes

    @staticmethod
    def order_psms(psms: List[PSM]) -> List[PSM]:
        psms.sort(key=lambda x: x.mz)
        return psms

    def add(self, psm: PSM) -> None:
        i = bisect(self.mz_list, psm.mz)
        self.tree.insert(i, psm)
        self.mz_list.insert(i, psm.mz)

    def update(self, psms: List[PSM]) -> None:
        for psm in psms:
            self.add(psm)

    def remove(self, psm: PSM):
        i = self.tree.index(psm)
        del (self.tree[i])
        del (self.mz_list[i])

    def _search(self, mz_boundary: Boundary, rt_boundary: Boundary, ook0_boundary: Boundary) -> List[PSM]:
        start_i = bisect_left(self.mz_list, mz_boundary.lower)
        if start_i == len(self.mz_list):
            return []
        start_mz = self.mz_list[start_i]
        if start_mz > mz_boundary.upper:
            return []
        res = []
        for i in range(start_i, len(self.mz_list)):
            psm = self.tree[i]
            if psm_attributes_in_bound(psm.mz, psm.rt, psm.ook0, mz_boundary, rt_boundary, ook0_boundary):
                res.append(psm)
            if self.mz_list[i] > mz_boundary.upper:
                break
        return res

    def get(self, mz: float, rt: float, ook0: float) -> List[PSM]:
        pass

    @property
    def psms(self) -> List[PSM]:
        return self.tree


@dataclass
class PsmSortedLinkedList(PsmSortedList):
    """
    Sorted Linked List performs exceptionally quickly in add time and, most importantly, search time.
    It is our first choice.
    Sorted Linked List has 2 properties: a tree composed of a list of PSM's,
    & a second list of mz values only which matches indices.
    """
    tree: List[PSM] = field(default_factory=lambda: deque())  #
    mz_list: [float] = field(default_factory=lambda: deque())  # this is a list of just mz values, for matching indexes


@dataclass
class PsmList(PsmTree):
    """
    List is costly to search but quick to add.
    """
    tree: List[PSM] = field(default_factory=lambda: [])

    @staticmethod
    def order_psms(psms: List[PSM]) -> List[PSM]:
        return psms

    def _search(self, mz_Boundary: Boundary, rt_Boundary: Boundary, ook0_Boundary: Boundary) -> List[PSM]:
        matches = []
        for psm in self.tree:
            if psm_attributes_in_bound(psm.mz, psm.rt, psm.ook0, mz_Boundary, rt_Boundary, ook0_Boundary):
                matches.append(psm)
        return matches

    def add(self, psm: PSM) -> None:
        self.tree.append(psm)

    def update(self, psms: List[PSM]) -> None:
        for psm in psms:
            self.add(psm)

    def remove(self, psm: PSM) -> None:
        self.tree.remove(psm)

    def get(self, mz: float, rt: float, ook0: float) -> List[PSM]:
        return [psm for psm in self.tree if psm.mz == mz and psm.rt == rt and psm.ook0 == ook0]

    @property
    def psms(self) -> List[PSM]:
        return self.tree

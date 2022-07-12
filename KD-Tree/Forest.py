from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from threading import Lock
from typing import Any, List, Dict
from bisect import bisect, insort, bisect_left
from collections import deque
from pykdtree import KDTree
from dataclasses import dataclass
from bintrees import FastBinaryTree, FastAVLTree, FastRBTree  # pip install bintrees
from .boundary import Boundary, psm_attributes_in_bound
from .point_util import *
from .red_black_tree import *

"""
PSM (becomes a python dict)
"{"mz":1000, "rt":250, "ook0":0.9, "Sequence": "PEPTIDE", "score": 2.0 ...}"
{"mz":1000, "rt":250, "ook0":0.9, "Sequence": "PEPTIDE", "score": 2.0 ...}
"""

"""
class PSM:
    def __init__(self, charge: int, mz: float, rt: float, ook0: float, sequence: str):
        self.charge = charge
        self.mz = mz
        self.rt = rt
        self.ook0 = ook0
        self.sequence = sequence
psm = PSM(charge=1, mz=100, rt=100, ook0=0.5, sequence="PEPTIDE")
"""

@dataclass
class PSM:
    charge: int
    mz: float
    rt: float
    ook0: float
    sequence: str

    def in_bounds(self, mz_bounds: Boundary, ook0_bounds: Boundary, rt_bounds: Boundary) -> bool:
        return mz_bounds.lower <= self.mz <= mz_bounds.upper and \
               ook0_bounds.lower <= self.ook0 <= ook0_bounds.upper and \
               rt_bounds.lower <= self.rt <= rt_bounds.upper

    def serialize(self) -> str:
        return f"{self.charge},{self.mz},{self.rt},{self.ook0},{self.sequence}\n"

    @staticmethod
    def deserialize(line: str) -> 'PSM':
        line_elems = line.rstrip().split(",")
        psm = PSM(int(line_elems[0]),
                  float(line_elems[1]),
                  float(line_elems[2]),
                  float(line_elems[3]),
                  str(line_elems[4]))
        return psm


#psm = PSM(charge=1, mz=100, rt=100, ook0=0.5, sequence="PEPTIDE")
""" 
    "charge": charge
    "mz": mz
    "rt": rt
    "ook0": ook0
    "Sequence": "PEPTIDE"
    "score": score
"""

"""
Trees:
Abstract (the structure & functions every tree should have)
KDTree
RedBlackTree
SortedList
List
"""
@dataclass
class AbstractPsmTree(ABC):
    tree: Any
    _lock: Lock = Lock()

    @abstractmethod
    def search(self, mz_bounds: Boundary, ook0_bounds: Boundary, rt_bounds: Boundary) -> List[PSM]:
        """
        searches the PSMTree over a given boundary. Return all psm's within the bounds
        """
        pass

    @abstractmethod
    def add(self, psm: PSM) -> None:
        """
        adds psm to interval tree with passed dimensions
        """
        pass

    @abstractmethod
    def len(self) -> int:
        pass

    @abstractmethod
    def save(self, file_name: str) -> None:
        """
        saves all psm's within tree to a text file
        """
        pass

    @abstractmethod
    def load(self, file_name: str) -> None:
        """
        adds psm's from text file to PSMTree
        """
        pass

"""
KDTree is our golden child. Thank you, KDTree.
"""
@dataclass
class PsmKDTree(AbstractPsmTree):
    tree: KDTree = KDTree()

    def search(self, mz_bounds: Boundary, ook0_bounds: Boundary, rt_bounds: Boundary) -> List[PSM]:
        lower_bounds = [mz_bounds.lower, ook0_bounds.lower, rt_bounds.lower]
        upper_bounds = [mz_bounds.upper, ook0_bounds.upper, tr_bounds.upper]
        self.tree.searchRange(self, List[PSM], lower_bounds, upper_bounds)

    def add(self, psm: PSM) -> None:
        pass

    def len(self) -> int:
        return len(self.tree)

    def save(self, file_name: str) -> None:
        with open(file_name, "w") as file:
            for psm in self.tree:
                file.write(psm.serialize())

    def load(self, file_name: str) -> None:
        with open(file_name) as file:
            for line in file:
                psm = PSM.deserialize(line)
                self.add(psm)


"""
RedBlackTree was crappy and thus abandoned
"""
@dataclass
class PsmRedBlackTree(AbstractPsmTree):
    tree: RedBlackTree = RedBlackTree()

    """
    search requests all important boundaries,
    utilizes red_black_tree.py's 'search_bounded' function (by sending only the upper & lower mz limits)
    which recursively collects & returns (to this 'search') a list of psms (called "data" in function).
    this list is then further filtered here (in 'search'), which returns a different list of best PSM matches.
    """
    def search(self, mz_bounds: Boundary, ook0_bounds: Boundary, rt_bounds: Boundary) -> List[PSM]:
        results = []
        for node in self.tree.search_bounded(mz_bounds.lower, mz_bounds.upper):
            for psm in node.data:
                if psm_attributes_in_bound(psm.mz, psm.ook0, psm.rt, mz_bounds, ook0_bounds, rt_bounds):
                    results.append(psm)
        return results

    def add(self, psm: PSM) -> None:
        node = self.tree.insert(psm.mz, [psm])
        if psm not in node.data:
            node.data.append(psm)

    def len(self) -> int:
        return len(self.tree)

    def save(self, file_name: str) -> None:
        with open(file_name, "w") as file:
            for psm in self.tree:
                file.write(psm.serialize())

    def load(self, file_name: str) -> None:
        with open(file_name) as file:
            for line in file:
                psm = PSM.deserialize(line)
                self.add(psm)


"""
PsmSortedList can use either binary search bisect functions or a deque to be much faster than PsmList but still long
bisect had significantly shorter search times but horrendous add times
deque had slightly better add times but horrendous search times (only after linked lists of 100k psm's)
"""
@dataclass
class PsmSortedList(AbstractPsmTree):
    tree: List[PSM] = field(default_factory=lambda:deque())
    mz_list:[float] = field(default_factory=lambda:deque())    #mz list for PSM

    def add(self, psm: PSM) -> None:
        i = bisect(self.mz_list, psm.mz)
        self.tree.insert(i, psm)
        self.mz_list.insert(i, psm.mz)

    def search(self, mz_bounds: Boundary, ook0_bounds: Boundary, rt_bounds: Boundary) -> List[PSM]:
        start_i = bisect_left(self.mz_list, mz_bounds.lower)
        if start_i == len(self.mz_list):
            return []
        #print(start_i, len(self.mz_list))
        start_mz = self.mz_list[start_i]
        if start_mz > mz_bounds.upper:
            return []
        res = []
        for i in range(start_i, len(self.mz_list)):
            psm = self.tree[i]
            if psm_attributes_in_bound(psm.mz, psm.ook0, psm.rt, mz_bounds, ook0_bounds, rt_bounds):
                res.append(psm)
            if self.mz_list[i] > mz_bounds.upper:
                break
        return res

    def len(self) -> int:
        return len(self.tree)

    def save(self, file_name: str) -> None:
        with open(file_name, "w") as file:
            for psm in self.tree:
                file.write(psm.serialize())

    def load(self, file_name: str) -> None:
        with open(file_name) as file:
            for line in file:
                psm = PSM.deserialize(line)
                self.add(psm)

"""
PsmList is a horrible choice. Do not choose this other than to see a worst case scenario. 
"""
@dataclass
class PsmList(AbstractPsmTree):
    tree: List[PSM] = field(default_factory=lambda:[])

    def search(self, mz_bounds: Boundary, ook0_bounds: Boundary, rt_bounds: Boundary) -> List[PSM]:
        matches = []
        for psm in self.tree:
            if psm_attributes_in_bound(psm.mz, psm.ook0, psm.rt, mz_bounds, ook0_bounds, rt_bounds):
                matches.append(psm)
        return matches


    """
    Adds a PSM to the "tree" list in a sorted manner.
    Utilizes binary search and tracks first, last, and self.mid indices. 
    Base Case: List is empty
    Case 1: New PSM is equal to first psm.
    Case 2: New PSM
    """
    def add(self, psm: PSM) -> None:
        if len(self.tree) == 0:
            self.tree.append(psm)
            return

    def len(self) -> int:
        return len(self.tree)

    def save(self, file_name: str) -> None:
        with open(file_name, "w") as file:
            for psm in self.tree:
                file.write(psm.serialize())

    def load(self, file_name: str) -> None:
        with open(file_name) as file:
            for line in file:
                psm = PSM.deserialize(line)
                self.add(psm)


@dataclass
class PsmBinaryTree(PsmAbstractTree):
    tree: FastBinaryTree = FastBinaryTree()

    def add(self, psm: PSM):
        self.tree.insert(psm.mz, psm)

    def search(self, psm, mz_bounds: Bounds, ook0_bounds: Bounds, rt_bounds: Bounds):
        res = self.tree[mz_bounds.lower:mz_bounds.upper].values()
        res = [psm for psm in res if psm.in_bounds(mz_bounds, ook0_bounds, rt_bounds)]
        return res

    def save(self, file_name: str) -> None:
        with open(file_name, "w") as file:
            for psm in self.tree:
                file.write(psm.serialize())

    def load(self, file_name: str) -> None:
        with open(file_name) as file:
            for line in file:
                psm = PSM.deserialize(line)
                self.add(psm)


@dataclass
class PsmAvlTree(PsmAbstractTree):
    tree: AVLTree = FastAVLTree()

    def add(self, psm: PSM):
        self.tree.insert(psm.mz, psm)

    def search(self, psm, mz_bounds: Bounds, ook0_bounds: Bounds, rt_bounds: Bounds):
        res = self.tree[mz_bounds.lower:mz_bounds.upper].values()
        res = [psm for psm in res if psm.in_bounds(mz_bounds, ook0_bounds, rt_bounds)]
        return res

    def save(self, file_name: str) -> None:
        with open(file_name, "w") as file:
            for psm in self.tree:
                file.write(psm.serialize())

    def load(self, file_name: str) -> None:
        with open(file_name) as file:
            for line in file:
                psm = PSM.deserialize(line)
                self.add(psm)

@dataclass
class PsmRBTree(PsmAbstractTree):
    tree: FastRBTree = FastRBTree()

    def add(self, psm: PSM):
        self.tree.insert(psm.mz, psm)

    def search(self, psm, mz_bounds: Bounds, ook0_bounds: Bounds, rt_bounds: Bounds):
        res = self.tree[mz_bounds.lower:mz_bounds.upper].values()
        res = [psm for psm in res if psm.in_bounds(mz_bounds, ook0_bounds, rt_bounds)]
        return res

    def save(self, file_name: str) -> None:
        with open(file_name, "w") as file:
            for psm in self.tree:
                file.write(psm.serialize())

    def load(self, file_name: str) -> None:
        with open(file_name) as file:
            for line in file:
                psm = PSM.deserialize(line)
                self.add(psm)
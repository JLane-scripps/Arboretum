from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from threading import Lock
from typing import Any, List, Dict
from bisect import bisect, insort, bisect_left
from collections import deque
from dataclasses import dataclass
from bintrees import FastBinaryTree, FastAVLTree, FastRBTree, BinaryTree, AVLTree, RBTree
import _pickle as cPickle
from .boundary import *
from .constants import *
from .point_util import *
from .red_black_tree import *
from .ty_kdtree import *

"""
A PSM is a Peptide Sequence Match. It's a set of values indicating charge-to-mass ratio (mz),
retention time, and a one-over-k-0 value. 
It also has a sequence, which is the literal peptide sequence, saved as an object and not a value used for search.
"""
@dataclass
class PSM:
    charge: int
    mz: float
    rt: float
    ook0: float
    sequence: str
    # psm = PSM(charge=1, mz=100, rt=100, ook0=0.5, sequence="PEPTIDE")

    def in_Boundary(self, mz_Boundary: Boundary, rt_Boundary: Boundary, ook0_Boundary: Boundary) -> bool:
        return mz_Boundary.lower <= self.mz <= mz_Boundary.upper and \
               rt_Boundary.lower <= self.rt <= rt_Boundary.upper and \
               ook0_Boundary.lower <= self.ook0 <= ook0_Boundary.upper

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
    def search(self, mz_Boundary: Boundary, rt_Boundary: Boundary, ook0_Boundary: Boundary) -> List[PSM]:
        """
        searches the PSMTree over a given boundary. Return all psm's within the Boundary
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

    """
    saves all psm's within tree to a text file
    """
    @abstractmethod
    def save(self, FILE_NAME):
        with open(FILE_NAME, "wb") as output_file:
            cPickle.dump(self.tree, output_file)

    """
    adds psm's from text file to PSMTree
    """
    @abstractmethod
    def load(self, FILE_NAME):
        with open(FILE_NAME, "rb") as input_file:
            self.tree = cPickle.load(input_file)

"""
KDTree, my beloved.
"""
@dataclass
class PsmKdTree(AbstractPsmTree):
    tree: KdTree = KdTree(3, [])

    def add(self, psm: PSM):
        self.tree.add(Point.create([psm.mz, psm.rt, psm.ook0], psm))

    def search(self, mz_bounds: Boundary, rt_bounds: Boundary, ook0_bounds: Boundary):
        bounds = [[mz_bounds.lower, mz_bounds.upper], [rt_bounds.lower, rt_bounds.upper], [ook0_bounds.lower, ook0_bounds.upper]]
        results = self.tree.get_bounded(bounds)
        return [res.data for res in results]

    def len(self) -> int:
        return len(self.tree)

    def save(self, file_name: str) -> None:
        with open(file_name, "w") as file:
            for point in self.tree._points:
                file.write(point.data.serialize())

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

    def search(self, mz_Boundary: Boundary, rt_Boundary: Boundary, ook0_Boundary: Boundary) -> List[PSM]:
        start_i = bisect_left(self.mz_list, mz_Boundary.lower)
        if start_i == len(self.mz_list):
            return []
        #print(start_i, len(self.mz_list))
        start_mz = self.mz_list[start_i]
        if start_mz > mz_Boundary.upper:
            return []
        res = []
        for i in range(start_i, len(self.mz_list)):
            psm = self.tree[i]
            if psm_attributes_in_bound(psm.mz, psm.rt, psm.ook0, mz_Boundary, ook0_Boundary, rt_Boundary):
                res.append(psm)
            if self.mz_list[i] > mz_Boundary.upper:
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

    def search(self, mz_Boundary: Boundary, rt_Boundary: Boundary, ook0_Boundary: Boundary) -> List[PSM]:
        matches = []
        for psm in self.tree:
            if psm_attributes_in_bound(psm.mz, psm.rt, psm.ook0, mz_Boundary, ook0_Boundary, rt_Boundary):
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
class AbstractPsmTree:
    tree: any

    def save(self, FILE_NAME):
        with open(FILE_NAME, "wb") as output_file:
            cPickle.dump(self.tree, output_file)

    def load(self, FILE_NAME):
        with open(FILE_NAME, "rb") as input_file:
            self.tree = cPickle.load(input_file)


@dataclass
class PsmBinTrees(AbstractPsmTree):

    def add(self, psm: PSM):
        self.tree.set_default(psm.mz, []).append(psm)

    def search(self, mz_bounds: Boundary, rt_bounds: Boundary, ook0_bounds: Boundary):
        res = self.tree[mz_bounds.lower:mz_bounds.upper+0.00001].values()
        return [psm for psm_list in res for psm in psm_list if psm.in_Boundary(mz_bounds, rt_bounds, ook0_bounds)]

    def clear(self):
        self.tree.clear()


@dataclass
class PsmFastBinaryTree(PsmBinTrees):
    tree: FastBinaryTree = FastBinaryTree()


@dataclass
class PsmFastAvlTree(PsmBinTrees):
    tree: FastAVLTree = FastAVLTree()


@dataclass
class PsmFastRBTree(PsmBinTrees):
    tree: FastRBTree = FastRBTree()


@dataclass
class PsmBinaryTree(PsmBinTrees):
    tree: BinaryTree = BinaryTree()


@dataclass
class PsmAvlTree(PsmBinTrees):
    tree: AVLTree = AVLTree()


@dataclass
class PsmRBTree(PsmBinTrees):
    tree: RBTree = RBTree()


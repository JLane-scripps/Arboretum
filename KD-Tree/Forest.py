from abc import ABC, abstractmethod
import numpy as np
from dataclasses import field
from threading import Lock
from typing import List
from bisect import bisect, bisect_left
from collections import deque
from intervaltree import IntervalTree
from bintrees import FastBinaryTree, FastAVLTree, FastRBTree, BinaryTree, AVLTree, RBTree
import _pickle as cPickle
from .boundary import *
from .point_util import *
from .red_black_tree import *
from .ty_kdtree import *




@dataclass
class PSM:
    """
    A PSM is a Peptide Sequence Match. It's a set of values indicating ion charge (charge),
    charge-to-mass ratio (mz), retention time (rt), and a one-over-k-0 value (ook0).
    It also has a string (sequence), the literal peptide sequence saved as an object and not a value used for search.
    It is ALWAYS listed in this order within this code for sake of consistency.
    """
    charge: int
    mz: float
    rt: float
    ook0: float
    sequence: str
    # Example: psm = PSM(charge=1, mz=100, rt=100, ook0=0.5, sequence="PEPTIDE")

    def in_boundary(self, mz_boundary: Boundary, rt_boundary: Boundary, ook0_boundary: Boundary) -> bool:
        """
        ensures each value is within acceptable range, regardless of how many dimensions are used for search.
        """
        return mz_boundary.lower <= self.mz <= mz_boundary.upper and \
               rt_boundary.lower <= self.rt <= rt_boundary.upper and \
               ook0_boundary.lower <= self.ook0 <= ook0_boundary.upper

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


@dataclass
class AbstractPsmTree(ABC):
    """
    The abstract idea of a tree that each of our tree types should follow, even when they are not true trees.
    """
    tree: Any
    _lock: Lock = Lock()

    @abstractmethod
    def search(self, mz_boundary: Boundary, rt_boundary: Boundary, ook0_boundary: Boundary) -> List[PSM]:
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

    def len(self) -> int:
        """
        returns the length of the tree
        """
        pass

    def save(self, FILE_NAME):
        """
        saves all psm's within tree to a text file
        """
        with open(FILE_NAME, "wb") as output_file:
            cPickle.dump(self.tree, output_file)

    def load(self, FILE_NAME):
        """
        adds psm's from text file to PSMTree
        """
        with open(FILE_NAME, "rb") as input_file:
            self.tree = cPickle.load(input_file)


@dataclass
class PsmKdTree(AbstractPsmTree):
    """
    KDTree, my beloved. Might need some thorough investigation and refinement. Currently slightly slower than SortedList
    """
    tree: KdTree = KdTree(3, [])

    def add(self, psm: PSM):
        self.tree.add(Point.create([psm.mz, psm.rt, psm.ook0], psm))

    def search(self, mz_bounds: Boundary, rt_bounds: Boundary, ook0_bounds: Boundary):
        bounds = [[mz_bounds.lower, mz_bounds.upper],
                  [rt_bounds.lower, rt_bounds.upper],
                  [ook0_bounds.lower, ook0_bounds.upper]]
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


@dataclass
class PsmIntervalTree(AbstractPsmTree):
    """
    Standard Interval Tree. Performs O(n * log n)
    """
    tree: IntervalTree = IntervalTree()
    ppm: int = 50

    def add(self, psm: PSM):
        ppm_offset = psm.mz * self.ppm / 1_000_000
        self.tree[psm.mz - ppm_offset:psm.mz + ppm_offset] = psm

    def search(self, mz_bounds: Boundary, rt_bounds: Boundary, ook0_bounds: Boundary):
        mz = np.mean([mz_bounds.lower, mz_bounds.upper])
        return [psm.data for psm in self.tree[mz] if psm.data.in_boundary(mz_bounds, rt_bounds, ook0_bounds)]

    def clear(self):
        self.tree.clear()

    def len(self) -> int:
        return len(self.tree)


@dataclass
class PsmSortedList(AbstractPsmTree):
    """
    PsmSortedList can use either binary search bisect functions or a deque to be much faster than PsmList...
    but still long. bisect had significantly shorter search times but horrendous add times
    deque had slightly better add times but horrendous search times (specifically after linked lists of 100k psm's)
    """
    tree: List[PSM] = field(default_factory=lambda: deque())  # The tree is composed of lists of PSM's
    mz_list: [float] = field(default_factory=lambda: deque())  # this is a list of just mz values, for matching indexes

    def add(self, psm: PSM) -> None:
        i = bisect(self.mz_list, psm.mz)
        self.tree.insert(i, psm)
        self.mz_list.insert(i, psm.mz)

    def search(self, mz_boundary: Boundary, rt_boundary: Boundary, ook0_boundary: Boundary) -> List[PSM]:
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
class PsmList(AbstractPsmTree):
    """
    PsmList is a horrible, messy choice. Do not choose this other than to see a worst case scenario.
    """
    tree: List[PSM] = field(default_factory=lambda: [])

    def search(self, mz_Boundary: Boundary, rt_Boundary: Boundary, ook0_Boundary: Boundary) -> List[PSM]:
        matches = []
        for psm in self.tree:
            if psm_attributes_in_bound(psm.mz, psm.rt, psm.ook0, mz_Boundary, rt_Boundary, ook0_Boundary):
                matches.append(psm)
        return matches

    def add(self, psm: PSM) -> None:
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
class PsmBinTrees(AbstractPsmTree):
    """
    Binary Search Tree. 1 Dimensional, so only compares mz values against each other for storage and search.
    Slightly above mediocre speed.
    Basis for some tree types, such as RB, AVL, Binary, and each of their "Fast" versions.
    """

    # pass in a PSM. The mz is extracted and put into list form, then appended to the tree.
    def add(self, psm: PSM):
        self.tree.set_default(psm.mz, []).append(psm)

    def search(self, mz_bounds: Boundary, rt_bounds: Boundary, ook0_bounds: Boundary):
        """
        pass in each of the 3 ranges. Only uses mz to sort in 1 dimension, but ensures all values are acceptable.
        hardcoded upper bound to prevent the upper bound itself from being sliced off and discounted.
        """
        res = self.tree[mz_bounds.lower:mz_bounds.upper + 0.00001].values()
        return [psm for psm_list in res for psm in psm_list if psm.in_boundary(mz_bounds, rt_bounds, ook0_bounds)]

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

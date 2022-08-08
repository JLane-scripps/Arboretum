from abc import ABC, abstractmethod
from enum import Enum
import numpy as np
from dataclasses import field, dataclass
from threading import Lock
from typing import List, Any
from bisect import bisect, bisect_left
from collections import deque
from kdtree import KdTree
from intervaltree import IntervalTree
import _pickle as pickle
from bintrees import AVLTree, RBTree, BinaryTree, FastAVLTree, FastBinaryTree, FastRBTree
from boundary import Boundary, psm_attributes_in_bound, get_mz_bounds, get_rt_bounds, get_ook0_bounds
from kdtree.point import Point
import sys

sys.setrecursionlimit(10 ** 6)


# All available tree types are as follows. This list must be updated whenever TreeTypes are removed or added.
class TreeType(Enum):
    KD_TREE = 1
    SORTED_LIST = 2
    LIST = 3
    FAST_AVL = 4
    FAST_RB = 5
    FAST_BINARY = 6
    BINARY = 7
    AVL = 8
    RB_TREE = 9
    INTERVAL_TREE = 10



# TreeType assignments corresponding to above
def psm_tree_constructor(tree_type: TreeType):
    if tree_type == TreeType.KD_TREE:
        return PsmKdTree()
    if tree_type == TreeType.SORTED_LIST:
        return PsmSortedList()
    if tree_type == TreeType.LIST:
        return PsmList()
    if tree_type == TreeType.BINARY:
        return PsmBinaryTree()
    if tree_type == TreeType.AVL:
        return PsmAvlTree()
    if tree_type == TreeType.RB_TREE:
        return PsmRBTree()
    if tree_type == TreeType.FAST_BINARY:
        return PsmFastBinaryTree()
    if tree_type == TreeType.FAST_AVL:
        return PsmFastAVLTree()
    if tree_type == TreeType.FAST_RB:
        return PsmFastRBTree()
    if tree_type == TreeType.INTERVAL_TREE:
        return PsmIntervalTree()
    else:
        return NotImplemented


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

    def __len__(self) -> int:
        """
        returns the length of the tree
        """
        return len(self.tree)

    def save(self, FILE_NAME):
        """
        saves all psm's within tree to a text file
        """
        with open(FILE_NAME, "wb") as output_file:
            pickle.dump(self.tree, FILE_NAME)

    def load(self, FILE_NAME):
        """
        adds psm's from text file to PSMTree
        """
        with open(FILE_NAME, "rb") as input_file:
            self.tree = pickle.load(FILE_NAME)


@dataclass
class PsmKdTree(AbstractPsmTree):
    """
    KDTree, my beloved. Slightly slower add time than IntervalTree (BEFORE 80k psm's ONLY). Better > 80K.
    Fastest Tree in the West.
    """
    tree: KdTree = field(default_factory=lambda: KdTree(3, []))

    def add(self, psm: PSM):
        p = Point([psm.mz, psm.rt, psm.ook0])
        p.data = psm
        self.tree.add_point(p)

    def search(self, mz_bounds: Boundary, rt_bounds: Boundary, ook0_bounds: Boundary):
        bounds = [[mz_bounds.lower, mz_bounds.upper],
                  [rt_bounds.lower, rt_bounds.upper],
                  [ook0_bounds.lower, ook0_bounds.upper]]
        results = self.tree.get_points_within_bounds(bounds)
        return [res.data for res in results]

    def __len__(self) -> int:
        return len(self.tree._points)


@dataclass
class PsmIntervalTree(AbstractPsmTree):
    """
    Standard Interval Tree. Performs O(n * log n).
    2nd best Tree, beating SortedList's disastrous add time.
    """
    tree: IntervalTree = field(default_factory=lambda: IntervalTree())
    ppm: int = 50

    def add(self, psm: PSM):
        ppm_offset = psm.mz * self.ppm / 1_000_000
        self.tree[psm.mz - ppm_offset:psm.mz + ppm_offset] = psm

    def search(self, mz_bounds: Boundary, rt_bounds: Boundary, ook0_bounds: Boundary):
        mz = np.mean([mz_bounds.lower, mz_bounds.upper])
        return [psm.data for psm in self.tree[mz] if psm.data.in_boundary(mz_bounds, rt_bounds, ook0_bounds)]

    def clear(self):
        self.tree.clear()

    def __len__(self) -> int:
        return len(self.tree)


@dataclass
class PsmSortedList(AbstractPsmTree):
    """
    Rank #3.
    PsmSortedList can use either binary search bisect functions or a deque to be much faster than PsmList...
    but still long. bisect had significantly shorter search times but horrendous add times
    deque had horrendous search times (after linked lists of 100k psm's) but slightly faster add times.
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
    Basis for tree types RB, AVL, & Binary
    """

    # pass in a PSM. The mz is extracted and put into list form, then appended to the tree.
    def add(self, psm: PSM):
        self.tree.set_default(psm.mz, []).append(psm)

    def search(self, mz_bounds: Boundary, rt_bounds: Boundary, ook0_bounds: Boundary):
        """
        pass in each of the 3 ranges. Only uses mz to sort in 1 dimension, but ensures all values are acceptable.
        hardcoded upper bound to prevent the upper bound itself from being sliced off and discounted.
        """
        res = self.tree.range_query_values([mz_bounds.lower, mz_bounds.upper])
        return [psm for psm_list in res for psm in psm_list if psm.in_boundary(mz_bounds, rt_bounds, ook0_bounds)]

    def clear(self):
        self.tree.clear()


    def __len__(self):
        return len(self.tree)


@dataclass
class PsmBinaryTree(PsmBinTrees):
    tree: BinaryTree = field(default_factory=lambda: BinaryTree())


@dataclass
class PsmAvlTree(PsmBinTrees):
    tree: AVLTree = field(default_factory=lambda: AVLTree())


@dataclass
class PsmRBTree(PsmBinTrees):
    tree: RBTree = field(default_factory=lambda: RBTree())


@dataclass
class PsmFastRBTree(PsmBinTrees):
    tree: FastRBTree = field(default_factory=lambda: FastRBTree())


@dataclass
class PsmFastAVLTree(PsmBinTrees):
    tree: FastAVLTree = field(default_factory=lambda: FastAVLTree())


@dataclass
class PsmFastBinaryTree(PsmBinTrees):
    tree: FastBinaryTree = field(default_factory=lambda: FastBinaryTree())

import ast
import sys
import numpy as np
from abc import ABC, abstractmethod
from bisect import bisect, bisect_left
from collections import deque
from dataclasses import field, dataclass
from enum import Enum
from intervaltree import IntervalTree
from kdtree import KdTree
from kdtree.point import Point
from typing import List, Any, Union

try:
    import cPickle as pickle
except:
    import pickle

from bintrees import BinaryTree, FastBinaryTree, AVLTree, FastAVLTree, RBTree, FastRBTree
from boundary import Boundary, psm_attributes_in_bound

sys.setrecursionlimit(10 ** 6)


# All available tree types are as follows. This list must be updated whenever TreeTypes are removed or added.
class TreeType(Enum):
    KD = 1
    SORTED_LIST = 2
    LIST = 3
    FAST_AVL = 4
    FAST_RB = 5
    FAST_BINARY = 6
    BINARY = 7
    AVL = 8
    RB = 9
    INTERVAL = 10
    SORTED_LINKED_LIST = 11


# TreeType assignments corresponding to above
def psm_tree_constructor(tree_type: Union[TreeType, str]):
    if tree_type == TreeType.KD or tree_type == 'kd_tree':
        return PsmKdTree()
    elif tree_type == TreeType.SORTED_LIST or tree_type == 'sorted_list':
        return PsmSortedList()
    elif tree_type == TreeType.LIST or tree_type == 'list':
        return PsmList()
    elif tree_type == TreeType.BINARY or tree_type == 'binary':
        return PsmBinaryTree()
    elif tree_type == TreeType.AVL or tree_type == 'avl':
        return PsmAvlTree()
    elif tree_type == TreeType.RB or tree_type == 'rb':
        return PsmRBTree()
    elif tree_type == TreeType.FAST_BINARY or tree_type == 'fast_binary':
        return PsmFastBinaryTree()
    elif tree_type == TreeType.FAST_AVL or tree_type == 'fast_avl':
        return PsmFastAVLTree()
    elif tree_type == TreeType.FAST_RB or tree_type == 'fast_rb':
        return PsmFastRBTree()
    elif tree_type == TreeType.INTERVAL or tree_type == 'interval':
        return PsmIntervalTree()
    elif tree_type == TreeType.SORTED_LINKED_LIST or tree_type == 'sorted_link_list':
        return PsmSortedLinkedList()
    else:
        raise Exception("Tree type not supported")


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
    data: dict

    # Example: psm = PSM(charge=1, mz=100, rt=100, ook0=0.5, sequence="PEPTIDE")

    def in_boundary(self, mz_boundary: Boundary, rt_boundary: Boundary, ook0_boundary: Boundary) -> bool:
        """
        ensures each value is within acceptable range, regardless of how many dimensions are used for search.
        """
        return mz_boundary.lower <= self.mz <= mz_boundary.upper and \
               rt_boundary.lower <= self.rt <= rt_boundary.upper and \
               ook0_boundary.lower <= self.ook0 <= ook0_boundary.upper

    def serialize(self) -> str:
        return f"{self.charge},{self.mz},{self.rt},{self.ook0},{self.data}\n"

    @staticmethod
    def deserialize(line: str) -> 'PSM':
        line_elems = line.rstrip().split(",")
        psm = PSM(int(line_elems[0]),
                  float(line_elems[1]),
                  float(line_elems[2]),
                  float(line_elems[3]),
                  ast.literal_eval(line_elems[4]))
        return psm


@dataclass
class AbstractPsmTree(ABC):
    """
    The abstract schematic of a tree that each of our tree types should follow, even when they are not true trees.
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
            pickle.dump(self.tree, output_file, -1)

    def load(self, FILE_NAME):
        """
        adds psm's from text file to PSMTree
        """
        with open(FILE_NAME, "rb") as input_file:
            self.tree = pickle.load(input_file)


@dataclass
class PsmKdTree(AbstractPsmTree):
    """

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
    PsmSortedList can use either binary search bisect functions or a deque to be much faster than PsmList...
    but still long. Bisect had significantly shorter search times but horrendous add times
    Deque had horrendous search times (after linked lists of 100k psm's) but slightly faster add times.
    """
    tree: List[PSM] = field(default_factory=lambda: list())  # The tree is composed of lists of PSM's
    mz_list: [float] = field(default_factory=lambda: list())  # this is a list of just mz values, for matching indexes

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

    def save(self, FILE_NAME):
        with open(FILE_NAME, "wb") as output_file:
            pickle.dump(self.tree, output_file)

    def load(self, FILE_NAME):
        with open(FILE_NAME, "rb") as input_file:
            self.tree = pickle.load(input_file)
        self.mz_list = [psm.mz for psm in self.tree]


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
class PsmList(AbstractPsmTree):
    """
    List is costly to search but quick to add.
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

    TODO: Update bintrees pickling process. Currently it calls add on every item. Get inspiration from interval_tree...
    """

    def add(self, psm: PSM):
        """
        Pass in a whole PSM. The mz is extracted and put into list form, then appended to the tree.
        """
        self.tree.set_default(psm.mz, []).append(psm)

    def search(self, mz_bounds: Boundary, rt_bounds: Boundary, ook0_bounds: Boundary):
        """
        pass in each of the 3 ranges as Boundary objects.
        Only uses mz to sort (because it's only sortable in 1 dimension), but ensures all values are acceptable.
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

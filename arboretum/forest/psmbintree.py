from dataclasses import dataclass, field
from typing import List

from bintrees.abctree import update_queue

from boundary import Boundary
from forest.psmtree import PsmTree
from psm import PSM
from bintrees import BinaryTree, FastBinaryTree, AVLTree, FastAVLTree, RBTree, FastRBTree


@dataclass
class PsmBinTree(PsmTree):
    """
    Binary Search Tree. 1 Dimensional, so only compares mz values against each other for storage and search.
    Slightly above mediocre speed.
    Basis for tree types RB, AVL, & Binary
    """

    @staticmethod
    def order_psms(psms: List[PSM]) -> List[PSM]:
        psms.sort(key=lambda x: x.mz)
        return update_queue(psms)

    def add(self, psm: PSM):
        self.tree.set_default(psm.mz, []).append(psm)

    def update(self, psms: List[PSM]) -> None:
        for psm in psms:
            self.add(psm)

    def remove(self, psm: PSM) -> None:
        psms = self.tree.get(psm.mz)

        if len(psms) == 1:
            self.tree.remove(psm.mz)
        else:
            psms.remove(psm)

    def _search(self, mz_bounds: Boundary, rt_bounds: Boundary, ook0_bounds: Boundary):
        """
        pass in each of the 3 ranges as Boundary objects.
        Only uses mz to sort (because it's only sortable in 1 dimension), but ensures all values are acceptable.
        """
        res = self.tree.range_query_values([mz_bounds.lower, mz_bounds.upper])
        return [psm for psm_list in res for psm in psm_list if psm.in_boundary(mz_bounds, rt_bounds, ook0_bounds)]

    def get(self, mz: float, rt: float, ook0: float) -> List[PSM]:
        psms = self.tree.get(mz)
        if psms is None:
            raise ValueError(f'no psm found with mz: {mz}')
        return [psm for psm in psms if psm.mz == mz and psm.rt == rt and psm.ook0 == ook0]

    @property
    def psms(self) -> List[PSM]:
        return [psm for res in self.tree.values() for psm in res]

    def clear(self):
        self.tree.clear()

    def __len__(self):
        return sum([len(l) for l in self.tree.values()])


@dataclass
class PsmBinaryTree(PsmBinTree):
    tree: BinaryTree = field(default_factory=lambda: BinaryTree())


@dataclass
class PsmAvlTree(PsmBinTree):
    tree: AVLTree = field(default_factory=lambda: AVLTree())


@dataclass
class PsmRBTree(PsmBinTree):
    tree: RBTree = field(default_factory=lambda: RBTree())


@dataclass
class PsmFastRBTree(PsmBinTree):
    tree: FastRBTree = field(default_factory=lambda: FastRBTree())


@dataclass
class PsmFastAVLTree(PsmBinTree):
    tree: FastAVLTree = field(default_factory=lambda: FastAVLTree())


@dataclass
class PsmFastBinaryTree(PsmBinTree):
    tree: FastBinaryTree = field(default_factory=lambda: FastBinaryTree())

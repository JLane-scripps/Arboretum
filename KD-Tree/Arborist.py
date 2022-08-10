import os
from dataclasses import dataclass, field
from threading import Lock
from typing import Dict
from Forest import TreeType, psm_tree_constructor, PSM, AbstractPsmTree, PsmKdTree, PsmIntervalTree, PsmList, \
    PsmSortedList, PsmSortedLinkedList, PsmBinTrees, PsmBinaryTree, PsmFastBinaryTree, PsmAvlTree, PsmFastAVLTree, \
    PsmRBTree, PsmFastRBTree
from boundary import Boundary
import shutil


@dataclass
class PSMArborist:
    """
    The Arborist is the handler of each tree. This class "tends to the forest".
    Stores PSM's in separate Trees (of TreeType) based on charge.
    Each Tree will have 1 consistent charge; therefore, expect 1 - 5 trees per TreeType.
    New Trees will be created for each new charge encountered
    """
    tree_type: TreeType
    trees: Dict[int, AbstractPsmTree] = field(default_factory=dict)

    _lock: Lock = Lock()

#  Find way to create a folder during runtime and save all trees within
    def save(self, directory):

        # Check if directory exists & remove if it does
        if os.path.exists(directory):
            shutil.rmtree(directory)

        print("Directory '% s' created" % directory)
        os.makedirs(directory)  # make new directory

        self._lock.acquire()
        for charge, tree in self.trees.items():
            file_name = f"{charge}.pkl"
            tree.save(os.path.join(directory, file_name))  # save tree's to charge.pkl
        self._lock.release()

    # pass a folder, look inside for saved files, and load them all as trees
    def load(self, directory):

        files = os.listdir(directory)

        self._lock.acquire()
        for file in files:
            charge = int(os.path.splitext(file)[0])  # /path/to/file.pkl -> file
            pms_tree = psm_tree_constructor(self.tree_type)
            pms_tree.load(os.path.join(directory, file))
            self.trees[charge] = pms_tree
        self._lock.release()

    def add(self, psm: PSM):
        """
        Adds a psm to the currently-used tree type, to the tree of correct charge.
        All trees less than 3 Dimensions prioritize sorting by mz limits.
        """
        charge = psm.charge
        if charge not in self.trees:
            self.trees[charge] = psm_tree_constructor(self.tree_type)
        self.trees[charge].add(psm)

    def search(self, charge, mz_bounds, rt_bounds, ook0_bounds):
        if charge not in self.trees:
            return []

        mz_bounds = Boundary(mz_bounds[0], mz_bounds[1])
        rt_bounds = Boundary(rt_bounds[0], rt_bounds[1])
        ook0_bounds = Boundary(ook0_bounds[0], ook0_bounds[1])
        results = self.trees[charge].search(mz_bounds, rt_bounds, ook0_bounds)

        return results

    def __len__(self):
        return sum([len(self.trees[charge]) for charge in self.trees])

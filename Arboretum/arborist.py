"""
-------------- Arborist --------------
A class designed to plant, grow, and tend
to "trees" of PSM's as they are received
by a mass spectrometer.
Supports a wide array of data storage
formats ("trees types").
--------------------------------------
"""

import os
from dataclasses import dataclass, field
from threading import Lock
from typing import Dict, Union
from forest import TreeType, psm_tree_constructor, PSM, AbstractPsmTree
from boundary import Boundary
import shutil


@dataclass
class PSMArborist:
    """
    The Arborist is the handler of each tree. This class "tends to the forest".
    Stores PSM's in separate Trees (of TreeType) based on charge.
    Each Tree will have 1 consistent charge; therefore, expect 1 - 5 trees per TreeType.
    New Trees will be created for each new charge encountered
    The chosen TreeType, declared here, can be specific or "TreeType" to allow all options and specified elsewhere.
    """
    tree_type: Union[TreeType, str] = TreeType.SORTED_LIST
    trees: Dict[int, AbstractPsmTree] = field(default_factory=dict)

    _lock: Lock = Lock()

    def save(self, directory):
        """
        Create a directory folder (name passed in) during runtime and save all trees within
        """
        # Check if directory exists & remove it if it does
        if os.path.exists(directory):
            shutil.rmtree(directory)

        # make new directory
        print("Directory '% s' created" % directory)
        os.makedirs(directory)

        self._lock.acquire()
        for charge, tree in self.trees.items():
            file_name = f"{charge}.pkl"
            tree.save(os.path.join(directory, file_name))  # save trees as [charge].pkl (i.e. "1.pkl", "2.pkl", etc)
        self._lock.release()


    def load(self, directory):
        """
        pass a folder, look inside for saved files, and load them all as trees
        """
        files = os.listdir(directory)

        self._lock.acquire()
        for file in files:
            charge = int(os.path.splitext(file)[0])  # /path/to/file.pkl -> file
            sapling = psm_tree_constructor(self.tree_type)
            sapling.load(os.path.join(directory, file))
            self.trees[charge] = sapling
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

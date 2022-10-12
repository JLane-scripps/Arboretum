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
import shutil

from forest import PsmTree, TreeType, psm_tree_constructor
from boundary import Boundary, get_mz_bounds, get_rt_bounds, get_ook0_bounds
from psm import PSM


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
    trees: Dict[int, PsmTree] = field(default_factory=dict)

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
            file_name = f"{charge}.txt"
            tree.save(os.path.join(directory, file_name),
                      as_pickle=False)  # save trees as [charge].pkl (i.e. "1.pkl", "2.pkl", etc)
        self._lock.release()

    def load(self, directory):
        """
        pass a folder, look inside for saved files, and load them all as trees
        """
        files = os.listdir(directory)

        self._lock.acquire()
        for file in files:
            charge = int(os.path.splitext(file)[0])  # /path/to/file.pkl -> file
            tree = psm_tree_constructor(self.tree_type)
            tree.load(os.path.join(directory, file), as_pickle=False)
            self.trees[charge] = tree
        self._lock.release()

    def add(self, charge: int, mz: float, rt: float, ook0: float, data: dict):
        """
        Adds a psm to the currently-used tree type, to the tree of correct charge.
        All trees less than 3 Dimensions prioritize sorting by mz limits.
        """
        psm = PSM(charge=charge, mz=mz, rt=rt, ook0=ook0, data=data)
        if psm.charge not in self.trees:
            self.trees[psm.charge] = psm_tree_constructor(self.tree_type)
        self.trees[psm.charge].add(psm)

    def search(self, charge: int, mz: float, rt: float, ook0: float, ppm: float, rt_offset: float,
               ook0_tolerance: float):
        if charge not in self.trees:
            return []

        mz_bounds = get_mz_bounds(mz, ppm)
        rt_bounds = get_rt_bounds(rt, rt_offset)
        ook0_bounds = get_ook0_bounds(ook0, ook0_tolerance)
        results = self.trees[charge]._search(mz_bounds, rt_bounds, ook0_bounds)
        return results

    def remove(self, charge: int, mz: float, rt: float, ook0: float, data: dict):
        psm = PSM(charge=charge, mz=mz, rt=rt, ook0=ook0, data=data)
        if psm.charge not in self.trees:
            raise ValueError(f'PSM not found. No tree with charge {charge}')
        self.trees[psm.charge].remove(psm)

    def __len__(self):
        return sum([len(self.trees[charge]) for charge in self.trees])

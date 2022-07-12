from dataclasses import dataclass, field
from enum import Enum
from threading import Lock
from typing import Dict, List
import _pickle as cPickle
from .Forest import *
from .boundary import *
from .constants import *
import random
import string


class TreeType(Enum):
    KD_TREE = 1
    SORTED_LIST = 2
    LIST = 3
    FAST_BINARY = 4
    FAST_AVL = 5
    FAST_RB = 6


def psm_tree_constructor(tree_type: TreeType):
    #if tree_type == TreeType.KD_TREE:
        #return PsmKDTree()
    if tree_type == TreeType.SORTED_LIST:
        return PsmSortedList()
    if tree_type == TreeType.LIST:
        return PsmList()
    else:
        return NotImplemented

#NEED TO ADD: Functions that calculate real Boundary values.
@dataclass
class PSMArborist:
    """
    Stores PSM's in separate Interval Trees according to charge
    New interval Trees will be created for each new charge encountered
    """
    tree_type: TreeType
    trees: Dict[int, AbstractPsmTree] = field(default_factory=dict)

    mz_ppm: int = 20
    ook0_tolerance: float = 0.05
    rt_offset: int = 100

    _latest_ms2_id: int = None
    _lock: Lock = Lock()

    def save(self, file):

        print(f"Saving PsmTree: {file}")
        self._lock.acquire()
        tmp_int_dict ={charge:self.trees[charge].tree for charge in self.trees}
        with open(file, "wb") as output_file:
            cPickle.dump(tmp_int_dict, output_file)
        self._lock.release()

    def load(self, file):
        print(f"loading PsmTree: {file}")
        self._lock.acquire()
        with open(file, "rb") as input_file:
            for charge, tree in cPickle.load(input_file).items():
                self.trees[charge] = AbstractPsmTree(tree)
        self._lock.release()

    def add(self, psm: Dict):
        """
        adds a psm dict to the interval tree according to the ~~interval: [mz-mz*ppm: mz+mz*ppm]~~  THIS IS OLD
        """
        charge = psm[PSM_CHARGE_KEY]
        if charge not in self.trees:
            self.trees[charge] = psm_tree_constructor(self.tree_type)
            self.trees[charge] = AbstractPsmTree(RedBlackTree())
        self.trees[charge].add(psm, ppm = self.mz_ppm)
        self._latest_ms2_id = psm['ms2_id']

    def search_psm(self, psm):  # mostly for testing
        return self.search(psm[PSM_CHARGE_KEY], psm[PSM_MZ_KEY], psm[PSM_OOK0_KEY], psm[PSM_RT_KEY])

    def search(self, charge, mz, ook0, rt, mz_ppm=None, ook0_tolerance=None, rt_offset=None) -> Dict:

        if charge not in self.trees:
            return []

        if mz_ppm is None:
            mz_ppm = self.mz_ppm

        if ook0_tolerance is None:
            ook0_tolerance = self.ook0_tolerance

        if rt_offset is None:
            rt_offset = self.rt_offset

        results = self.trees[charge].search(mz, ook0, rt, mz_ppm, ook0_tolerance, rt_offset)

        return results

    def len(self):
        return sum([self.trees[charge].len() for charge in self.trees])
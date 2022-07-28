from enum import Enum
from typing import Dict
from .boundary import Boundary
from .constants import *
from .Forest import*



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

    def save(self, file):
        print(f"Saving PsmTree: {file}")
        self._lock.acquire()
        with open(file, "wb") as output_file:
            cPickle.dump(self.trees, output_file)
        self._lock.release()

    def load(self, file):
        print(f"loading PsmTree: {file}")
        self._lock.acquire()
        with open(file, "rb") as input_file:
            self.trees = cPickle.load(input_file)
        self._lock.release()

    def add(self, psm: PSM):
        """
        Adds a psm to the currently-used tree type, to the tree of correct charge.
        All trees less than 3 Dimensions prioritize sorting by mz limits.
        """
        charge = psm.charge
        if charge not in self.trees:
            self.trees[charge] = psm_tree_constructor(self.tree_type)()
        self.trees[charge].add(psm)

    def search(self, charge, mz_bounds, rt_bounds, ook0_bounds):
        if charge not in self.trees:
            return []

        results = self.trees[charge].search(mz_bounds, rt_bounds, ook0_bounds)

        return results

    def __len__(self):
        return sum([len(self.trees[charge]) for charge in self.trees])

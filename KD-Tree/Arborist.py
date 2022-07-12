from enum import Enum
import _pickle as cPickle
from .Forest import*


class TreeType(Enum):
    KD_TREE = 1
    SORTED_LIST = 2
    LIST = 3
    FAST_BINARY = 4
    FAST_AVL = 5
    FAST_RB = 6
    BINARY = 7
    AVL = 8
    RB_TREE = 9


def psm_tree_constructor(tree_type: TreeType):
    if tree_type == TreeType.KD_TREE:
        return PsmKdTree
    if tree_type == TreeType.SORTED_LIST:
        return PsmSortedList
    if tree_type == TreeType.LIST:
        return PsmList
    if tree_type == TreeType.FAST_BINARY:
        return PsmFastBinaryTree
    if tree_type == TreeType.FAST_AVL:
        return PsmFastAvlTree
    if tree_type == TreeType.FAST_RB:
        return PsmFastRBTree
    if tree_type == TreeType.BINARY:
        return PsmBinaryTree
    if tree_type == TreeType.AVL:
        return PsmAvlTree
    if tree_type == TreeType.RB_TREE:
        return PsmRBTree
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
    rt_offset: int = 100
    ook0_tolerance: float = 0.05

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
            self.trees[charge] = AbstractPsmTree(self.tree_type)
        self.trees[charge].add(psm, ppm = self.mz_ppm)
        self._latest_ms2_id = psm['ms2_id']

    def search_psm(self, psm):  # mostly for testing
        return self.search(psm[PSM_CHARGE_KEY], psm[PSM_MZ_KEY], psm[PSM_RT_KEY], psm[PSM_OOK0_KEY])

    def search(self, charge, mz, rt, ook0, mz_ppm=None, rt_offset=None, ook0_tolerance=None) -> Dict:

        if charge not in self.trees:
            return []

        if mz_ppm is None:
            mz_ppm = self.mz_ppm

        if rt_offset is None:
            rt_offset = self.rt_offset

        if ook0_tolerance is None:
            ook0_tolerance = self.ook0_tolerance

        results = self.trees[charge].search(mz, rt, ook0, mz_ppm, rt_offset, ook0_tolerance)

        return results

    def len(self):
        return sum([self.trees[charge].len() for charge in self.trees])
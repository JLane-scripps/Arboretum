import math
from dataclasses import dataclass, field
from typing import List

from boundary import Boundary
from forest.psmtree import PsmTree
from psm import PSM


def convert_to_int(mz:float, precision:int, floor=True):
    if floor == True:
        return math.floor(round(mz, precision)*(10**precision))
    elif floor == False:
        return math.ceil(round(mz, precision)*(10**precision))


@dataclass
class PsmHashtable(PsmTree):
    """
    Standard Interval Tree. Performs O(n * log n).
    """
    tree: dict = field(default_factory=lambda: {})
    precision: int = 2


    @staticmethod
    def order_psms(psms: List[PSM]) -> List[PSM]:
        return psms

    def add(self, psm: PSM) -> None:
        key = convert_to_int(psm.mz, self.precision)
        self.tree.setdefault(key, []).append(psm)

    def remove(self, psm: PSM) -> None:
        key = convert_to_int(psm.mz, self.precision)
        self.tree[key].remove(psm)

    def _search(self, mz_boundary: Boundary, rt_boundary: Boundary, ook0_boundary: Boundary) -> List[PSM]:

        lower_key = convert_to_int(mz_boundary.lower, self.precision)
        upper_key = convert_to_int(mz_boundary.upper, self.precision, floor=False)

        psms = []
        for key in range(lower_key, upper_key+1):
            if key in self.tree:
                psms.extend(self.tree[key])

        return [psm for psm in psms if psm.in_boundary(mz_boundary, rt_boundary, ook0_boundary)]

    def get(self, mz: float, rt: float, ook0: float) -> List[PSM]:
        pass

    @property
    def psms(self) -> List[PSM]:
        return list(self.tree)
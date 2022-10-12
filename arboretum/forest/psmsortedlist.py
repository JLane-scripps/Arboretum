import math
from dataclasses import dataclass, field
from typing import List

from sortedcontainers import SortedDict

from boundary import Boundary
from forest.psmtree import PsmTree
from psm import PSM


@dataclass
class PsmSortedList(PsmTree):
    tree: SortedDict = field(default_factory=lambda: SortedDict())

    @staticmethod
    def order_psms(psms: List[PSM]) -> List[PSM]:
        psms.sort(key=lambda x: x.mz)
        return psms

    def add(self, psm: PSM) -> None:
        if psm.mz not in self.tree:
            self.tree[psm.mz] = []
        self.tree[psm.mz].append(psm)

    def remove(self, psm: PSM) -> None:
        if psm.mz not in self.tree:
            raise ValueError
        self.tree[psm.mz].remove(psm)

    def search(self, mz_boundary: Boundary, rt_boundary: Boundary, ook0_boundary: Boundary) -> List[PSM]:
        keys = self.tree.irange(mz_boundary.lower, mz_boundary.upper)
        psms = []
        for key in keys:
            psms.extend(self.tree[key])
        return [psm for psm in psms if psm.in_boundary(mz_boundary, rt_boundary, ook0_boundary)]

    def get(self, mz: float, rt: float, ook0: float) -> List[PSM]:
        if mz not in self.tree:
            raise ValueError

        return [psm for psm in self.tree[mz] if psm.mz == mz and psm.rt == rt and psm.ook0 == ook0]

    @property
    def psms(self) -> List[PSM]:
        return [psm for key in self.tree for psm in self.tree[key]]

    def __len__(self):
        return sum([len(self.tree[key]) for key in self.tree])

    def clear(self):
        self.tree.clear()


def convert_to_int(mz:float, precision:int, floor=True):
    if floor == True:
        return math.floor(round(mz, precision)*(10**precision))
    elif floor == False:
        return math.ceil(round(mz, precision)*(10**precision))

@dataclass
class PsmHashtable(PsmTree):
    tree: SortedDict = field(default_factory=lambda: SortedDict())
    precision: int = 2

    @staticmethod
    def order_psms(psms: List[PSM]) -> List[PSM]:
        psms.sort(key=lambda x: x.mz)
        return psms

    def add(self, psm: PSM) -> None:
        key = convert_to_int(psm.mz, self.precision)
        if key not in self.tree:
            self.tree[key] = []
        self.tree[key].append(psm)

    def remove(self, psm: PSM) -> None:
        key = convert_to_int(psm.mz, self.precision)
        if key not in self.tree:
            raise ValueError
        self.tree[key].remove(psm)

    def search(self, mz_boundary: Boundary, rt_boundary: Boundary, ook0_boundary: Boundary) -> List[PSM]:
        lower_key = convert_to_int(mz_boundary.lower, self.precision)
        upper_key = convert_to_int(mz_boundary.upper, self.precision, floor=False)

        psms = []
        for key in range(lower_key, upper_key+1):
            if key in self.tree:
                psms.extend(self.tree[key])

        return [psm for psm in psms if psm.in_boundary(mz_boundary, rt_boundary, ook0_boundary)]

    def get(self, mz: float, rt: float, ook0: float) -> List[PSM]:
        key = convert_to_int(mz, self.precision)
        if key not in self.tree:
            raise ValueError

        return [psm for psm in self.tree[key] if psm.mz == mz and psm.rt == rt and psm.ook0 == ook0]

    @property
    def psms(self) -> List[PSM]:
        return [psm for key in self.tree for psm in self.tree[key]]

    def __len__(self):
        return sum([len(self.tree[key]) for key in self.tree])

    def clear(self):
        self.tree.clear()


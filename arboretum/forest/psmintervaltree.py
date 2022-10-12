from dataclasses import dataclass, field
from typing import List

from intervaltree import IntervalTree

from boundary import Boundary
from forest.psmtree import PsmTree
from psm import PSM


@dataclass
class PsmIntervalTree(PsmTree):
    """
    Standard Interval Tree. Performs O(n * log n).
    """
    tree: IntervalTree = field(default_factory=lambda: IntervalTree())
    ppm: int = 50

    @staticmethod
    def order_psms(psms: List[PSM]) -> List[PSM]:
        return psms

    def add(self, psm: PSM):
        ppm_offset = psm.mz * self.ppm / 1_000_000
        self.tree[psm.mz - ppm_offset:psm.mz + ppm_offset] = psm

    def remove(self, psm: PSM):
        psms = [psm.data for psm in self.tree[psm.mz]]
        psms.remove(psm)

    def _search(self, mz_bounds: Boundary, rt_bounds: Boundary, ook0_bounds: Boundary):
        psms = [interval.data for interval in self.tree[mz_bounds.lower:mz_bounds.upper+0.0001]]  # make inclusive
        return [psm for psm in psms if psm.in_boundary(mz_bounds, rt_bounds, ook0_bounds)]

    def get(self, mz: float, rt: float, ook0: float) -> List[PSM]:
        psms = [interval.data for interval in self.tree[mz]]
        if len(psms) > 0 and psms[0].mz == mz and psms[0].rt == rt and psms[0].ook0 == ook0:
            return psms
        else:
            raise ValueError

    @property
    def psms(self) -> List[PSM]:
        return [interval.data for interval in self.tree.items()]

    def clear(self):
        self.tree.clear()

    def __len__(self) -> int:
        return len(self.tree)

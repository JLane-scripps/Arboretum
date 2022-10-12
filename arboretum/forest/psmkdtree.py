from dataclasses import dataclass, field
from typing import List

from kdtree import KdTree

from boundary import Boundary
from forest.psmtree import PsmTree
from point_util import Point
from psm import PSM

# TODO: Fix kd-tree package to support deletion

@dataclass
class PsmKdTree(PsmTree):
    """
    KD-Tree implementation of AbstractPsmTree
    """
    tree: KdTree = field(default_factory=lambda: KdTree(3, []))

    @staticmethod
    def order_psms(psms: List[PSM]) -> List[PSM]:
        return psms

    def add(self, psm: PSM):
        p = Point([psm.mz, psm.rt, psm.ook0])
        p.data = psm
        self.tree.add_point(p)

    def remove(self, psm: PSM):
        psm = self.get(psm.mz, psm.rt, psm.ook0)


    def search(self, mz_bounds: Boundary, rt_bounds: Boundary, ook0_bounds: Boundary):
        bounds = [[mz_bounds.lower, mz_bounds.upper],
                  [rt_bounds.lower, rt_bounds.upper],
                  [ook0_bounds.lower, ook0_bounds.upper]]
        results = self.tree.get_points_within_bounds(bounds)
        return [res.data for res in results]

    def get(self, mz: float, rt: float, ook0: float) -> List[PSM]:
        p = Point([mz, rt, ook0])
        psm = self.tree.get_nearest(p).data

        if psm.mz == mz and psm.rt == rt and psm.ook0 == ook0:
            return psm
        else:
            raise ValueError

    @property
    def psms(self) -> List[PSM]:
        return [point.data for point in self.tree._points]

    def __len__(self) -> int:
        return len(self.tree._points)

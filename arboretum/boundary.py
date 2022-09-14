from dataclasses import dataclass


@dataclass
class Boundary:
    __slots__ = "lower", "upper"
    lower: float
    upper: float


def psm_attributes_in_bound(mz, rt, ook0, mz_bounds: Boundary, rt_bounds: Boundary, ook0_bounds: Boundary) -> bool:
    return mz_bounds.lower <= mz <= mz_bounds.upper and \
           rt_bounds.lower <= rt <= rt_bounds.upper and \
           ook0_bounds.lower <= ook0 <= ook0_bounds.upper

"""
receives mass-to-charge ration (mz) as a float, and ppm as a float
uses mz as focal point, subtracting & adding (mz * ppm)
returns upper bound & lower bound
" -> Boundary " ensures it returns as class Boundary, from boundary.py
"""
def get_mz_bounds(mz: float, ppm: float) -> Boundary:
    return Boundary(lower=mz - mz * ppm / 1_000_000,
                    upper=mz + mz * ppm / 1_000_000)

"""
receives retention time (rt) as a float, & offset as a float
simple rt +- offset
returns lower bound & upper bound as a Boundary object
"""
def get_rt_bounds(rt: float, offset: float) -> Boundary:
    return Boundary(lower=rt - offset,
                    upper=rt + offset)

"""
receives one over k0 (ook0) as a float, and tolerance as a float
uses ook0 as focal point, for +- of (ook0 * tolerance)
returns lower bound & upper bound as a class Boundary object 
"""
def get_ook0_bounds(ook0: float, tolerance: float) -> Boundary:
    return Boundary(lower=ook0 - ook0 * tolerance,
                    upper=ook0 + ook0 * tolerance)

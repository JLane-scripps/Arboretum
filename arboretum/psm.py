import ast
from dataclasses import dataclass

from boundary import Boundary


@dataclass
class PSM:
    """
    A PSM is a Peptide Sequence Match. It's a set of values indicating ion charge (charge),
    charge-to-mass ratio (mz), retention time (rt), and a one-over-k-0 value (ook0).
    It also has a string (sequence), the literal peptide sequence saved as an object and not a value used for search.
    It is ALWAYS listed in this order within this code for sake of consistency.
    """
    charge: int
    mz: float
    rt: float
    ook0: float
    data: dict

    # Example: psm = PSM(charge=1, mz=100, rt=100, ook0=0.5, sequence="PEPTIDE")

    def in_boundary(self, mz_boundary: Boundary, rt_boundary: Boundary, ook0_boundary: Boundary) -> bool:
        """
        ensures each value is within acceptable range, regardless of how many dimensions are used for search.
        """
        return mz_boundary.lower <= self.mz <= mz_boundary.upper and \
               rt_boundary.lower <= self.rt <= rt_boundary.upper and \
               ook0_boundary.lower <= self.ook0 <= ook0_boundary.upper

    def serialize(self) -> str:
        return f"{self.charge},{self.mz},{self.rt},{self.ook0},{self.data}\n"

    @staticmethod
    def deserialize(line: str) -> 'PSM':
        line_elems = line.rstrip().split(",")
        psm = PSM(charge=int(line_elems[0]),
                  mz=float(line_elems[1]),
                  rt=float(line_elems[2]),
                  ook0=float(line_elems[3]),
                  data=ast.literal_eval(line_elems[4]))
        return psm

    def __eq__(self, other):
        return self.charge == other.charge and self.mz == other.mz and self.rt == other.rt and self.ook0 == other.ook0 \
               and self.data == other.data
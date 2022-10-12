from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, Union, List

from boundary import Boundary, get_mz_bounds, get_rt_bounds, get_ook0_bounds
from psm import PSM

try:
    import cPickle as pickle
except:
    import pickle


@dataclass
class PsmTree(ABC):
    """
    The abstract schematic of a tree that each of our tree types should follow, even when they are not true trees.
    """
    tree: Any
    __psms: Union[List[PSM], None] = None

    def __post_init__(self):
        if self.__psms:
            psms = self.order_psms(self.__psms)
            for psm in psms:
                self.tree.add(psm)
            del self.__psms

    @staticmethod
    @abstractmethod
    def order_psms(psms: List[PSM]) -> List[PSM]:
        """
        given a list of psms reorder them so so that tree insertion will not cause excess re-balancing
        """
        pass

    def tsearch(self, mz: float, rt: float, ook0: float, ppm: float, rt_offset: float, ook0_tolerance: float) -> List[PSM]:
        """
        tsearch : tolerance search. searches the tree with given tolerances.
        """
        mz_bounds = get_mz_bounds(mz, ppm)
        rt_bounds = get_rt_bounds(rt, rt_offset)
        ook0_bounds = get_ook0_bounds(ook0, ook0_tolerance)
        return self._search(mz_bounds, rt_bounds, ook0_bounds)

    def search(self, mz_boundary: List[float], rt_boundary: List[float], ook0_boundary: List[float]) -> List[PSM]:

        if len(mz_boundary) != 2 or len(rt_boundary) != 2 or len(ook0_boundary) != 2:
            raise ValueError('Incorrect boundary arguments. Boundary should contain two items: [lower, upper]')

        return self._search(Boundary(mz_boundary[0], mz_boundary[1]),
                     Boundary(rt_boundary[0], rt_boundary[1]),
                     Boundary(ook0_boundary[0], ook0_boundary[1]))

    @abstractmethod
    def _search(self, mz_boundary: Boundary, rt_boundary: Boundary, ook0_boundary: Boundary) -> List[PSM]:
        """
        searches the PSMTree over a given boundary. Return all psm's within the Boundary
        """
        pass

    @abstractmethod
    def get(self, mz: float, rt: float, ook0: float) -> List[PSM]:
        """
        Return all psm's with matching values
        """
        pass

    @abstractmethod
    def add(self, psm: PSM) -> None:
        """
        adds psm to tree with passed dimensions
        """
        pass

    def update(self, psms: List[PSM]) -> None:
        """
        adds all psm's to the tree
        """
        for psm in psms:
            self.add(psm)

    @abstractmethod
    def remove(self, psm: PSM) -> None:
        """
        removes psm from tree with passed dimensions. Throws ValueError if psm not in tree
        """
        pass

    def discard(self, psm: PSM) -> bool:
        """
        removes psm from tree with passed dimensions. Returns true if removed, false if not present
        """
        try:
            self.remove(psm)
        except ValueError:
            return False
        return True

    @property
    @abstractmethod
    def psms(self) -> List[PSM]:
        """
        returns psms in the tree
        """
        pass

    def __len__(self) -> int:
        """
        returns the length of the tree
        """
        return len(self.tree)

    def save(self, file_name: str, as_pickle: bool = False):
        """
        loads psms to file: either pkl or txt
        """
        if as_pickle:
            self.to_pickle(file_name)
        else:
            self.to_file(file_name)

    def load(self, file_name: str, as_pickle: bool = False):
        """
        loads psms from file: either pkl or txt
        """
        if as_pickle:
            self.from_pickle(file_name)
        else:
            self.from_file(file_name)

    def to_pickle(self, file_name: str):
        with open(file_name, "wb") as file:
            pickle.dump(self.tree, file, -1)

    def from_pickle(self, file_name: str):
        with open(file_name, "rb") as file:
            self.tree = pickle.load(file)

    def to_file(self, file_name: str):
        with open(file_name, "w") as file:
            for psm in self.psms:
                file.write(psm.serialize())

    def from_file(self, file_name: str):
        psms = []
        with open(file_name, "r") as file:
            for line in file:
                psms.append(PSM.deserialize(line))
        psms = self.order_psms(psms)
        self.update(psms)

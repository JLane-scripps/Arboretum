# All available tree types are as follows. This list must be updated whenever TreeTypes are removed or added.
from enum import Enum, auto


class TreeType(Enum):
    KD = auto()
    INTERVAL = auto()
    SORTED_LIST = auto()
    HASHTABLE = auto()
    HASHTABLE_MED = auto()
    HASHTABLE_LARGE = auto()
    FAST_AVL = auto()
    FAST_RB = auto()
    FAST_BINARY = auto()
    BINARY = auto()
    AVL = auto()
    RB = auto()

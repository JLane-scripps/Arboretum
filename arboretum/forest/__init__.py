import sys
from typing import Union

from forest.psmtree import PsmTree
from forest.psmbintree import PsmBinaryTree, PsmAvlTree, PsmRBTree, PsmFastBinaryTree, PsmFastAVLTree, PsmFastRBTree
from forest.psmintervaltree import PsmIntervalTree
from forest.psmkdtree import PsmKdTree
from forest.psmsortedlist import PsmSortedList, PsmHashtable
from forest.treetypes import TreeType

sys.setrecursionlimit(10 ** 6)


# TreeType assignments corresponding to above
def psm_tree_constructor(tree_type: Union[TreeType, str]):
    if tree_type == TreeType.KD or tree_type == 'kd_tree':
        # return PsmKdTree()
        return NotImplementedError
    elif tree_type == TreeType.BINARY or tree_type == 'binary':
        return PsmBinaryTree()
    elif tree_type == TreeType.AVL or tree_type == 'avl':
        return PsmAvlTree()
    elif tree_type == TreeType.RB or tree_type == 'rb':
        return PsmRBTree()
    elif tree_type == TreeType.FAST_BINARY or tree_type == 'fast_binary':
        return PsmFastBinaryTree()
    elif tree_type == TreeType.FAST_AVL or tree_type == 'fast_avl':
        return PsmFastAVLTree()
    elif tree_type == TreeType.FAST_RB or tree_type == 'fast_rb':
        return PsmFastRBTree()
    elif tree_type == TreeType.INTERVAL or tree_type == 'interval':
        # return PsmIntervalTree()
        return NotImplementedError
    elif tree_type == TreeType.SORTED_LIST or tree_type == 'sorted_list':
        return PsmSortedList()
    elif tree_type == TreeType.HASHTABLE or tree_type == 'hashtable':
        return PsmHashtable()
    elif tree_type == TreeType.HASHTABLE_MED or tree_type == 'hashtable_med':
        return PsmHashtable(precision=3)
    elif tree_type == TreeType.HASHTABLE_LARGE or tree_type == 'hashtable_large':
        return PsmHashtable(precision=4)
    else:
        raise Exception("Tree type not supported")
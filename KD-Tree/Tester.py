import unittest
import time
import random
from .Arborist import *


def generate_random_psm() -> PSM:
    letters = 'ARNDCEQGHILKMFPSTWYV'
    peptide_string = ''.join(random.choice(letters) for i in range(random.randint(6, 30)))
    return PSM(
        charge=random.randint(1, 3),
        mz=random.uniform(400, 1800),
        rt=random.uniform(0, 10_000),
        ook0=random.uniform(0.4, 1.8),
        sequence=peptide_string
    )


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
    INTERVAL_TREE = 10


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
    if tree_type == TreeType.INTERVAL_TREE:
        return PsmIntervalTree
    else:
        return NotImplemented

def test_by_psm_tree_type(tree_type: TreeType):
    class PsmTreeTester(unittest.TestCase):
        PPM = 20
        RT_OFF = 100
        OOK0_TOL = 0.05

        def setUp(self):
            self.setup_tree = psm_tree_constructor(tree_type)
            self.tree = self.setup_tree()
            self.psms = [PSM(1, 1005.0, 250, 0.9, "PEPTIDE"),
                         PSM(1, 1100.0, 260, 0.9, "PETIDE"),
                         PSM(1, 1150.0, 260, 0.9, "PETIDE"),
                         PSM(1, 1200.0, 252, 0.9, "PEP"),
                         PSM(2, 1000.0, 250, 0.9, "PEP"),
                         PSM(2, 1050.0, 300, 0.9, "PEPTI"),
                         PSM(3, 1000.0, 250, 0.9, "PEPTIDE"),
                         PSM(1, 3000.0, 250, 0.9, "PEPTIDE")
                         ]

        """
        If mass is the same, tree is not adding results
        each node can store a list of psm's
        (problem: now really intermingled with PSMTree)
        or 
        make a left-child added as the new psm.
        """

        def test_add_search_performance(self):
            print("running test_add_search_performance for", self.setup_tree)
            performance_dict = {}
            for n in [10, 100, 500, 1_000]:
                tree = self.setup_tree()
                start_time = time.time()
                psms = [generate_random_psm() for _ in range(n)]
                _ = [tree.add(psm) for psm in psms]
                add_time = (time.time() - start_time)
                print(f"N: {n:,}, add_time: {add_time:,}")

                start_time = time.time()
                _ = [self.assertTrue(len(tree.search(get_mz_bounds(psm.mz, PsmTreeTester.PPM),
                                                     get_rt_bounds(psm.rt, PsmTreeTester.RT_OFF),
                                                     get_ook0_bounds(psm.ook0, PsmTreeTester.OOK0_TOL))) > 0) for psm in psms]
                search_time = (time.time() - start_time)
                print(f"N: {n:,}, search_time: {search_time}")
                print(f"N: {n:,}, search_time_per_psm: {(time.time() - start_time) / n}")

                performance_dict[n] = {'add_time': add_time, 'search_time': search_time}
            print(performance_dict)

        def test_add(self):
            self.tree.add(self.psms[0])

        def test_search_basic(self):
            psm = self.psms[0]
            self.tree.add(psm)
            results = self.tree.search(get_mz_bounds(psm.mz, PsmTreeTester.PPM),
                                       get_rt_bounds(psm.rt, PsmTreeTester.RT_OFF),
                                       get_ook0_bounds(psm.ook0, PsmTreeTester.OOK0_TOL))
            self.assertTrue(psm in results)

        def test_search_all(self):
            for psm in self.psms:
                self.tree.add(psm)
            for psm in self.psms:
                results = self.tree.search(get_mz_bounds(psm.mz, PsmTreeTester.PPM),
                                           get_rt_bounds(psm.rt, PsmTreeTester.RT_OFF),
                                           get_ook0_bounds(psm.ook0, PsmTreeTester.OOK0_TOL))
                self.assertTrue(psm in results)
            # print(self.tree)

        def test_edges(self):
            for x in self.psms:
                self.tree.add(x)
            for psm in self.psms:
                results = self.tree.search(get_mz_bounds(psm.mz, 0),
                                           get_rt_bounds(psm.rt, 0),
                                           get_ook0_bounds(psm.ook0, 0))

                #if psm in results:
                #    print("psm success: ", psm)
                #elif psm not in results:
                #    print("psm not in results")

        def test_not_within(self):
            for psm in self.psms:
                self.tree.add(psm)
            for psm in self.psms:
                results = self.tree.search(get_mz_bounds(0.1, 0),
                                           get_rt_bounds(0.01, 0),
                                           get_ook0_bounds(0.001, 0))
                self.assertFalse(psm in results)

        """
        setting the lower bound of the range equal to the exact mz value 
        """

        def test_mz_lower(self):
            for new in self.psms:
                self.tree.add(new)
            for psm in self.psms:
                ppm_offset = psm.mz * 20 / 1_000_000
                mz_bounds = Boundary(psm.mz, psm.mz + ppm_offset)

                results = self.tree.search(mz_bounds,
                                           get_rt_bounds(psm.rt, 10000),
                                           get_ook0_bounds(psm.ook0, 1000))
                #print(results)
                self.assertTrue(psm in results)

        def test_mz_upper(self):
            for psm in self.psms:
                self.tree.add(psm)
            for psm in self.psms:
                ppm_offset = psm.mz * 20 / 1_000_000
                mz_bounds = Boundary(psm.mz - ppm_offset, psm.mz)

                results = self.tree.search(mz_bounds,
                                           get_rt_bounds(psm.rt, 10000),
                                           get_ook0_bounds(psm.ook0, 1000))
                #print(results)
                self.assertTrue(psm in results)

        def test_rt_lower(self):
            for psm in self.psms:
                self.tree.add(psm)
            for psm in self.psms:
                rt_offset = psm.rt + PsmTreeTester.RT_OFF
                rt_bounds = Boundary(psm.rt, psm.rt + rt_offset)

                results = self.tree.search(get_mz_bounds(psm.mz, 1000),
                                           rt_bounds,
                                           get_ook0_bounds(psm.ook0, 1000))
                self.assertTrue(psm in results)

        def test_rt_upper(self):
            for psm in self.psms:
                self.tree.add(psm)
            for psm in self.psms:
                rt_offset = psm.rt + PsmTreeTester.RT_OFF
                rt_bounds = Boundary(psm.rt - rt_offset, psm.rt)

                results = self.tree.search(get_mz_bounds(psm.mz, 1000),
                                           rt_bounds,
                                           get_ook0_bounds(psm.ook0, 1000))
                self.assertTrue(psm in results)

        def test_ook0_lower(self):
            for psm in self.psms:
                self.tree.add(psm)
            for psm in self.psms:
                ook0_tolerance = psm.ook0 + psm.ook0 * PsmTreeTester.OOK0_TOL
                ook0_bounds = Boundary(psm.ook0, psm.ook0 + ook0_tolerance)

                results = self.tree.search(get_mz_bounds(psm.mz, 1000),
                                           get_rt_bounds(psm.rt, 10000),
                                           ook0_bounds)
                self.assertTrue(psm in results)

        def test_ook0_upper(self):
            for psm in self.psms:
                self.tree.add(psm)
            for psm in self.psms:
                ook0_tolerance = psm.ook0 - psm.ook0 * PsmTreeTester.OOK0_TOL
                ook0_bounds = Boundary(psm.ook0 - ook0_tolerance, psm.ook0)

                results = self.tree.search(get_mz_bounds(psm.mz, 1000),
                                           get_rt_bounds(psm.rt, 10000),
                                           ook0_bounds)
                self.assertTrue(psm in results)

        def test_save_load(self):
            for psm in self.psms:
                self.tree.add(psm)
            self.tree.save('temp.txt')
            tree2 = self.setup_tree()
            tree2.load('temp.txt')
            for psm in self.psms:
                results = tree2.search(get_mz_bounds(psm.mz, PsmTreeTester.PPM),
                                       get_rt_bounds(psm.rt, PsmTreeTester.RT_OFF),
                                       get_ook0_bounds(psm.ook0, PsmTreeTester.OOK0_TOL))
                self.assertTrue(psm in results)

        # 1) add psms to tree
        # 2) search all psms's
        #   - verify each psm in results
        #   - verify that all psms returned by results are within bounds
        """def test_generate_random_psm(self):
            N = 10_000
            psms = [generate_random_psm() for i in range(N)]
            for psm in psms:
                self.tree.add(psm)
                results = self.tree.search(get_mz_bounds(psm.mz, PsmTreeTester.PPM),
                                           get_rt_bounds(psm.rt, PsmTreeTester.RT_OFF), 
                                           get_ook0_bounds(psm.ook0, PsmTreeTester.OOK0_TOL))
                self.assertTrue(psm in results)"""


    return PsmTreeTester


"""for tree_type in TreeType:
    class TreeTester(test_by_psm_tree_type(tree_type)):
        pass"""

class KDTreeTester(test_by_psm_tree_type(TreeType.KD_TREE)):
    pass

class SLTreeTester(test_by_psm_tree_type(TreeType.SORTED_LIST)):
    pass

class ListTreeTester(test_by_psm_tree_type(TreeType.LIST)):
    pass

class FastBinTreeTester(test_by_psm_tree_type(TreeType.FAST_BINARY)):
    pass

class FastAVLTreeTester(test_by_psm_tree_type(TreeType.FAST_AVL)):
    pass

class FastRBTreeTester(test_by_psm_tree_type(TreeType.FAST_RB)):
    pass

class BinTreeTester(test_by_psm_tree_type(TreeType.BINARY)):
    pass

class AVLTreeTester(test_by_psm_tree_type(TreeType.AVL)):
    pass

class RBTreeTester(test_by_psm_tree_type(TreeType.RB_TREE)):
    pass

class IntervalTreeTester(test_by_psm_tree_type(TreeType.INTERVAL_TREE)):
    pass

if __name__ == '__main__':
    unittest.main()

"""psms = []
        for i in range(10_000):
            psms.append(generate_random_psm())"""

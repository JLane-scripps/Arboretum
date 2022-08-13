import os
import unittest
import random
import time
from forest import TreeType, psm_tree_constructor, PSM, AbstractPsmTree, PsmKdTree, PsmIntervalTree, PsmList, \
    PsmSortedList, PsmSortedLinkedList, PsmBinTrees, PsmBinaryTree, PsmFastBinaryTree, PsmAvlTree, PsmFastAVLTree, \
    PsmRBTree, PsmFastRBTree
from boundary import Boundary, get_mz_bounds, get_rt_bounds, get_ook0_bounds


def generate_random_psm() -> PSM:
    letters = 'ARNDCEQGHILKMFPSTWYV'
    peptide_string = ''.join(random.choice(letters) for i in range(random.randint(6, 30)))
    return PSM(
        charge=random.randint(1, 5),
        mz=random.uniform(100, 1800),
        rt=random.uniform(0, 10_000),
        ook0=random.uniform(0.4, 1.8),
        sequence=peptide_string
    )


def test_by_psm_tree_type(tree_type: TreeType):
    class PsmTreeTester(unittest.TestCase):
        PPM = 50
        RT_OFF = 250
        OOK0_TOL = 0.05

        def setUp(self):
            self.tree = psm_tree_constructor(tree_type)
            self.psms = [PSM(1, 1005.0, 250, 0.9, "PEPTIDE"),
                         PSM(1, 1100.0, 260, 0.9, "PETIDE"),
                         PSM(1, 1150.0, 260, 0.9, "PETIDE"),
                         PSM(1, 1200.0, 252, 0.9, "PEP"),
                         PSM(2, 2000.0, 250, 0.9, "PEP"),
                         PSM(2, 1050.0, 300, 0.9, "PEPTI"),
                         PSM(3, 5000.0, 250, 0.9, "PEPTIDE"),
                         PSM(1, 3000.0, 250, 0.9, "PEPTIDE")
                         ]  # static psm's for the sake of testing.

        def test_add_search_performance(self):
            """
            # Utilizes every function.
            # Creates tree.
            # Adds static psms and searches for them.
            # Generates random psms in increasingly large quantities and searches for them.
            # Measures search times.
            # Saves & Loads trees.
            """
            print("\nPerformance for", tree_type.name)
            performance_dict = {}
            for n in [10, 100, 500, 1_000, 2_500, 5_000, 10_000, 50_000, 75_000, 100_000 , 250_000, 400_000, 500_000] : #]:
                tree = psm_tree_constructor(tree_type)
                start_time = time.time()
                psms = [generate_random_psm() for _ in range(n)]
                _ = [tree.add(psm) for psm in psms]
                add_time = (time.time() - start_time)
                print(f"N: {n:,}, add_time: {add_time:,}")

                start_time = time.time()
                _ = [self.assertTrue(len(tree.search(get_mz_bounds(psm.mz, PsmTreeTester.PPM),
                                                     get_rt_bounds(psm.rt, PsmTreeTester.RT_OFF),
                                                     get_ook0_bounds(psm.ook0, PsmTreeTester.OOK0_TOL))) > 0) for psm in
                     psms]
                search_time = (time.time() - start_time)
                print(f"N: {n:,}, search_time: {search_time}")
                print(f"N: {n:,}, search_time_per_psm: {(time.time() - start_time) / n}")

                performance_dict[n] = {'add_time': add_time, 'search_time': search_time}
            print(performance_dict)

            file_name = os.path.join('C:/Users/jeffl/OneDrive/Scripps Research Institute/Arboretum/Performance Logs '
                                     'Temp', tree_type.name +'.txt')
            with open(file_name, 'w') as newfile:
                newfile.write(str(performance_dict))
            newfile.close()


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

        def test_edges(self):
            for x in self.psms:
                self.tree.add(x)
            for psm in self.psms:
                results = self.tree.search(get_mz_bounds(psm.mz, 0),
                                           get_rt_bounds(psm.rt, 0),
                                           get_ook0_bounds(psm.ook0, 0))
                self.assertTrue(psm in results)

        def test_not_within(self):
            for psm in self.psms:
                self.tree.add(psm)
            for psm in self.psms:
                results = self.tree.search(get_mz_bounds(0.1, 0),
                                           get_rt_bounds(0.01, 0),
                                           get_ook0_bounds(0.001, 0))
                self.assertFalse(psm in results)


        """ ----------- Searching for each value's upper & lower bound ------------ """

        def test_mz_lower(self):
            for new in self.psms:
                self.tree.add(new)
            for psm in self.psms:
                ppm_offset = psm.mz * 20 / 1_000_000
                mz_bounds = Boundary(psm.mz, psm.mz + ppm_offset)
                results = self.tree.search(mz_bounds,
                                           get_rt_bounds(psm.rt, 10000),
                                           get_ook0_bounds(psm.ook0, 1000))
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

#ADD TIME STAMPS FOR SAVE & LOAD
        def test_save_load(self):
            for psm in self.psms:
                self.tree.add(psm)
            start_time = time.time()
            self.tree.save('temp.txt')
            save_time = time.time() - start_time
            print("save time: ", save_time)
            tree2 = psm_tree_constructor(tree_type)

            start_time = time.time()
            tree2.load('temp.txt')
            load_time = time.time() - start_time
            print ("load time: ", load_time)
            for psm in self.psms:
                results = tree2.search(get_mz_bounds(psm.mz, PsmTreeTester.PPM),
                                       get_rt_bounds(psm.rt, PsmTreeTester.RT_OFF),
                                       get_ook0_bounds(psm.ook0, PsmTreeTester.OOK0_TOL))
                self.assertTrue(psm in results)

    return PsmTreeTester


""" ----------- Repeats each function above for each respective Tree Type -------------- """

"""
class KDTreeTester(test_by_psm_tree_type(TreeType.KD_TREE)):
    pass

class SLTreeTester(test_by_psm_tree_type(TreeType.SORTED_LIST)):
    pass

class ListTreeTester(test_by_psm_tree_type(TreeType.LIST)):
    pass
    
class IntervalTreeTester(test_by_psm_tree_type(TreeType.INTERVAL_TREE)):
    pass
"""

class SLTreeTester(test_by_psm_tree_type(TreeType.SORTED_LIST)):
    pass

class IntervalTreeTester(test_by_psm_tree_type(TreeType.INTERVAL_TREE)):
    pass

"""
class BinTreeTester(test_by_psm_tree_type(TreeType.BINARY)):
    pass


class AVLTreeTester(test_by_psm_tree_type(TreeType.AVL)):
    pass


class RBTreeTester(test_by_psm_tree_type(TreeType.RB_TREE)):
    pass


class FastAVLTree(test_by_psm_tree_type(TreeType.FAST_AVL)):
    pass

class FastRBTree(test_by_psm_tree_type(TreeType.FAST_RB)):
    pass

class FastBinaryTree(test_by_psm_tree_type(TreeType.FAST_BINARY)):
    pass
"""

if __name__ == '__main__':
    unittest.main()

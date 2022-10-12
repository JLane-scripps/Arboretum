import os
import unittest
import random
import time

import numpy as np

from arboretum.forest import TreeType, psm_tree_constructor
from arboretum.boundary import Boundary, get_mz_bounds, get_rt_bounds, get_ook0_bounds
from psm import PSM


def generate_random_psm() -> PSM:
    letters = 'ARNDCEQGHILKMFPSTWYV'
    peptide_string = ''.join(random.choice(letters) for i in range(random.randint(6, 30)))
    mz = np.random.normal(1000, 10)
    return PSM(
        charge=random.randint(1, 5),
        mz=mz,
        rt=random.uniform(0, 250),
        ook0=mz/1000 + random.uniform(-0.2, 0.2),
        data={'sequence':peptide_string}
    )


def test_by_psm_tree_type(tree_type: TreeType):
    class PsmTreeTester(unittest.TestCase):
        PPM = 50
        RT_OFF = 100
        OOK0_TOL = 0.05

        def setUp(self):
            self.tree = psm_tree_constructor(tree_type)
            self.psms = [PSM(1, 1005.0, 250, 0.9, {'sequence':'PEPTIDE'}),
                         PSM(1, 1100.0, 260, 0.9, {'sequence':'PETIDE'}),
                         PSM(1, 1150.0, 260, 0.9, {'sequence':'PETIDE'}),
                         PSM(1, 1200.0, 252, 0.9, {'sequence':'PEP'}),
                         PSM(2, 2000.0, 250, 0.9, {'sequence':'PEP'}),
                         PSM(2, 1050.0, 300, 0.9, {'sequence':'PEPTI'}),
                         PSM(3, 5000.0, 250, 0.9, {'sequence':'PEPTIDE'}),
                         PSM(1, 3000.0, 250, 0.9, {'sequence':'PEPTIDE'})
                         ]  # static psm's for the sake of testing.

        def test_add(self):
            self.tree.add(self.psms[0])

        def test_add_negative(self):
            psm = PSM(-1, -1005.0, -250, -0.9, {'sequence':'PEPTIDE'})
            self.tree.add(psm)
            results = self.tree.search(get_mz_bounds(psm.mz, PsmTreeTester.PPM),
                                       get_rt_bounds(psm.rt, PsmTreeTester.RT_OFF),
                                       get_ook0_bounds(psm.ook0, PsmTreeTester.OOK0_TOL))
            self.assertTrue(psm in results)

        def test_add_zero(self):
            psm = PSM(0, 0,0,0, {'sequence':'PEPTIDE'})
            self.tree.add(psm)
            results = self.tree.search(get_mz_bounds(psm.mz, PsmTreeTester.PPM),
                                       get_rt_bounds(psm.rt, PsmTreeTester.RT_OFF),
                                       get_ook0_bounds(psm.ook0, PsmTreeTester.OOK0_TOL))
            self.assertTrue(psm in results)

        def test_remove(self):
            self.tree.add(self.psms[0])
            self.assertEqual(1, len(self.tree))
            self.tree.remove(self.psms[0])
            self.assertEqual(0, len(self.tree))

        def test_remove_all(self):
            for psm in self.psms:
                self.tree.add(psm)

            for i, psm in enumerate(self.psms, 1):
                self.tree.remove(psm)
                self.assertEqual(len(self.psms) - i, len(self.tree))

        def test_get(self):
            self.tree.add(self.psms[0])
            self.assertEqual(self.psms[0], self.tree.get(self.psms[0].mz, self.psms[0].rt, self.psms[0].ook0)[0])

        def test_get_fail(self):
            self.tree.add(self.psms[0])
            self.assertRaises(ValueError, self.tree.get, self.psms[1].mz, self.psms[1].rt, self.psms[1].ook0)

        def test_get_all(self):
            for psm in self.psms:
                self.tree.add(psm)

            for psm in self.psms:
                self.assertEqual(psm, self.tree.get(psm.mz, psm.rt, psm.ook0)[0])

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

        def test_dup(self):
            self.tree.add(self.psms[0])
            self.tree.add(self.psms[0])
            self.tree.add(self.psms[0])
            self.tree.add(self.psms[0])
            self.assertEqual(4, len(self.tree))

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

        def test_random(self):
            psms = [generate_random_psm() for i in range(1000)]
            for psm in psms:
                self.tree.add(psm)

            for psm in psms:
                results = self.tree.search(get_mz_bounds(psm.mz, PsmTreeTester.PPM),
                                           get_rt_bounds(psm.rt, PsmTreeTester.RT_OFF),
                                           get_ook0_bounds(psm.ook0, PsmTreeTester.OOK0_TOL))
                self.assertTrue(psm in results)


        # ADD TIME STAMPS FOR SAVE & LOAD
        def test_save_load(self):
            for psm in self.psms:
                self.tree.add(psm)
            start_time = time.time()
            self.tree.save('temp.txt')
            save_time = time.time() - start_time
            tree2 = psm_tree_constructor(tree_type)

            start_time = time.time()
            tree2.load('temp.txt')
            load_time = time.time() - start_time
            for psm in self.psms:
                results = tree2.search(get_mz_bounds(psm.mz, PsmTreeTester.PPM),
                                       get_rt_bounds(psm.rt, PsmTreeTester.RT_OFF),
                                       get_ook0_bounds(psm.ook0, PsmTreeTester.OOK0_TOL))
                self.assertTrue(psm in results)

    return PsmTreeTester


""" ----------- Repeats each function above for each respective Tree Type -------------- """


#class KDTreeTester(test_by_psm_tree_type(TreeType.KD)): pass

#class SLTreeTester(test_by_psm_tree_type(TreeType.SORTED_LIST)): pass

#class SLLTreeTester(test_by_psm_tree_type(TreeType.SORTED_LINKED_LIST)): pass

#class IntervalTreeTester(test_by_psm_tree_type(TreeType.INTERVAL)):pass

class PsmSortedList(test_by_psm_tree_type(TreeType.SORTED_LIST)):pass

class PsmHashtable(test_by_psm_tree_type(TreeType.HASHTABLE)):pass

class BinTreeTester(test_by_psm_tree_type(TreeType.BINARY)):pass

class AVLTreeTester(test_by_psm_tree_type(TreeType.AVL)):pass

class RBTreeTester(test_by_psm_tree_type(TreeType.RB)):pass

class FastAVLTree(test_by_psm_tree_type(TreeType.FAST_AVL)):pass

class FastRBTree(test_by_psm_tree_type(TreeType.FAST_RB)):pass

class FastBinaryTree(test_by_psm_tree_type(TreeType.FAST_BINARY)):pass

if __name__ == '__main__':
    unittest.main()

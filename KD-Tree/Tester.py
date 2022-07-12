import unittest
import time
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

def psm_tree_constructor(tree_type: TreeType):
    if tree_type == TreeType.KD_TREE:
        return PsmKDTree
    if tree_type == TreeType.SORTED_LIST:
        return PsmSortedList
    if tree_type == TreeType.LIST:
        return PsmList
    else:
        return NotImplemented

class PsmTreeTester(unittest.TestCase):
    PPM = 20
    OOK0_TOL = 0.05
    RT_OFF = 100


    def setUp(self):
        self.setup_tree = psm_tree_constructor(TreeType.KD_TREE)
        self.tree = self.setup_tree()
        self.psms = [PSM(1, 1005, 250, 0.9, "PEPTIDE"),
                     PSM(1, 1100, 260, 0.9, "PETIDE"),
                     PSM(1, 1150, 260, 0.9, "PETIDE"),
                     PSM(1, 1200, 252, 0.9, "PEP"),
                     PSM(2, 1000, 250, 0.9, "PEP"),
                     PSM(2, 1050, 300, 0.9, "PEPTI"),
                     PSM(3, 1000, 250, 0.9, "PEPTIDE"),
                     PSM(1, 3000, 250, 0.9, "PEPTIDE")
                     ]

    """
    If mass is the same, tree is not adding results
    each node can store a list of psm's
    (problem: now really intermignled with PSMTree)
    or 
    make a left-child added as the new psm.
    """

    def test_add(self):
        self.tree.add(self.psms[0])

    def test_search_basic(self):
        psm = self.psms[0]
        self.tree.add(psm)
        results = self.tree.search(get_mz_bounds(psm.mz, PsmTreeTester.PPM),
                                   get_ook0_bounds(psm.ook0, PsmTreeTester.OOK0_TOL),
                                   get_rt_bounds(psm.rt, PsmTreeTester.RT_OFF))
        self.assertTrue(psm in results)

    def test_search_all(self):
        for psm in self.psms:
            self.tree.add(psm)
        for psm in self.psms:
            results = self.tree.search(get_mz_bounds(psm.mz, PsmTreeTester.PPM),
                                       get_ook0_bounds(psm.ook0, PsmTreeTester.OOK0_TOL),
                                       get_rt_bounds(psm.rt, PsmTreeTester.RT_OFF))
            self.assertTrue(psm in results)
        # print(self.tree)

    def test_edges(self):
        for x in self.psms:
            self.tree.add(x)
        for psm in self.psms:
            results = self.tree.search(get_mz_bounds(psm.mz, 0),
                                       get_ook0_bounds(psm.ook0, 0),
                                       get_rt_bounds(psm.rt, 0))

            if psm in results: print("psm success: ", psm)
            elif psm not in results: print("psm not in results")

    def test_not_within(self):
        for psm in self.psms:
            self.tree.add(psm)
        for psm in self.psms:
            results = self.tree.search(get_mz_bounds(0.1, 0),
                                       get_ook0_bounds(0.001, 0),
                                       get_rt_bounds(0.01, 0))
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
                                       get_ook0_bounds(psm.ook0, 1000),
                                       get_rt_bounds(psm.rt, 10000))
            print(results)
            self.assertTrue(psm in results)

    def test_mz_upper(self):
        for psm in self.psms:
            self.tree.add(psm)
        for psm in self.psms:
            ppm_offset = psm.mz * 20 / 1_000_000
            mz_bounds = Boundary(psm.mz - ppm_offset, psm.mz)

            results = self.tree.search(mz_bounds,
                                       get_ook0_bounds(psm.ook0, 1000),
                                       get_rt_bounds(psm.rt, 10000))
            print(results)
            self.assertTrue(psm in results)

    def test_ook0_lower(self):
        for psm in self.psms:
            self.tree.add(psm)
        for psm in self.psms:
            ook0_tolerance = psm.ook0 + psm.ook0 * PsmTreeTester.OOK0_TOL
            ook0_bounds = Boundary(psm.ook0, psm.ook0 + ook0_tolerance)

            results = self.tree.search(get_mz_bounds(psm.mz, 1000),
                                       ook0_bounds,
                                       get_rt_bounds(psm.rt, 10000))
            self.assertTrue(psm in results)

    def test_ook0_upper(self):
        for psm in self.psms:
            self.tree.add(psm)
        for psm in self.psms:
            ook0_tolerance = psm.ook0 - psm.ook0 * PsmTreeTester.OOK0_TOL
            ook0_bounds = Boundary(psm.ook0 - ook0_tolerance, psm.ook0)

            results = self.tree.search(get_mz_bounds(psm.mz, 1000),
                                       ook0_bounds,
                                       get_rt_bounds(psm.rt, 10000))
            self.assertTrue(psm in results)

    def test_rt_lower(self):
        for psm in self.psms:
            self.tree.add(psm)
        for psm in self.psms:
            rt_offset = psm.rt + PsmTreeTester.RT_OFF
            rt_bounds = Boundary(psm.rt, psm.rt + rt_offset)

            results = self.tree.search(get_mz_bounds(psm.mz, 1000),
                                       get_ook0_bounds(psm.ook0, 1000),
                                       rt_bounds)
            self.assertTrue(psm in results)

    def test_rt_upper(self):
        for psm in self.psms:
            self.tree.add(psm)
        for psm in self.psms:
            rt_offset = psm.rt + PsmTreeTester.RT_OFF
            rt_bounds = Boundary(psm.rt - rt_offset, psm.rt)

            results = self.tree.search(get_mz_bounds(psm.mz, 1000),
                                       get_ook0_bounds(psm.ook0, 1000),
                                       rt_bounds)
            self.assertTrue(psm in results)

    def test_save_load(self):
        for psm in self.psms:
            self.tree.add(psm)
        self.tree.save('temp.txt')
        tree2 = PsmSortedList()
        tree2.load('temp.txt')
        for psm in self.psms:
            results = tree2.search(get_mz_bounds(psm.mz, PsmTreeTester.PPM),
                                   get_ook0_bounds(psm.ook0, PsmTreeTester.OOK0_TOL),
                                   get_rt_bounds(psm.rt, PsmTreeTester.RT_OFF))
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
                                       get_ook0_bounds(psm.ook0, PsmTreeTester.OOK0_TOL),
                                       get_rt_bounds(psm.rt, PsmTreeTester.RT_OFF))
            self.assertTrue(psm in results)"""

    def test_add_search_performance(self):
        print("running test_add_search_performance")
        performance_dict = {}
        for n in [1_000, 5_000, 10_000, 50_000, 100_000, 500_000, 1_000_000, 5_000_000]:
            tree = self.setup_tree()
            start_time = time.time()
            psms = [generate_random_psm() for _ in range(n)]
            _ = [tree.add(psm) for psm in psms]
            add_time = (time.time() - start_time)
            print(f"N: {n:,}, add_time: {add_time:,}")

            start_time = time.time()
            _ = [self.assertTrue(len(tree.search(get_mz_bounds(psm.mz, PsmTreeTester.PPM),
                                                 get_ook0_bounds(psm.ook0, PsmTreeTester.OOK0_TOL),
                                                 get_rt_bounds(psm.rt, PsmTreeTester.RT_OFF))) > 0) for psm in psms]
            search_time = (time.time() - start_time)
            print(f"N: {n:,}, search_time: {search_time}")
            print(f"N: {n:,}, search_time_per_psm: {(time.time() - start_time) / n}")

            performance_dict[n] = {'add_time': add_time, 'search_time': search_time}
        print(performance_dict)

    def test_search_performance2(self):
        print("running test_add_search_performance2")
        return
        performance_dict = {}
        tree = self.setup_tree()
        for n in range(100):
            psms = [generate_random_psm() for _ in range(100_000)]

            start_time = time.time()
            _ = [tree.add(psm) for psm in psms]
            add_time = (time.time() - start_time)
            print(f"N: {n}, add_time: {add_time}")

            start_time = time.time()
            _ = [self.assertTrue(len(tree.search(get_mz_bounds(psm.mz, PsmTreeTester.PPM),
                                                 get_ook0_bounds(psm.ook0, PsmTreeTester.OOK0_TOL),
                                                 get_rt_bounds(psm.rt, PsmTreeTester.RT_OFF))) > 0) for psm in psms]
            search_time = (time.time() - start_time)
            print(f"N: {n}, search_time: {search_time}")
            performance_dict[n] = {'add_time': add_time, 'search_time': search_time}
            print(performance_dict)


if __name__ == '__main__':
    unittest.main()

"""psms = []
        for i in range(10_000):
            psms.append(generate_random_psm())"""

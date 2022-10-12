import random
import time

import numpy as np
from matplotlib import pyplot as plt

from boundary import get_mz_bounds, get_rt_bounds, get_ook0_bounds
from forest import TreeType, psm_tree_constructor
from psm import PSM

num_points = 10
num_psms = 10_000
n = 1000

AMINOACIDS = 'ARNDCEQGHILKMFPSTWYV'

def generate_random_psm() -> PSM:
    letters = 'ARNDCEQGHILKMFPSTWYV'
    peptide_string = ''.join(random.choice(letters) for i in range(random.randint(6, 30)))
    mz = np.random.normal(1000, 250)
    return PSM(
        charge=random.randint(1, 5),
        mz=mz,
        rt=random.uniform(0, 5000),
        ook0=mz/1000 + random.uniform(-0.2, 0.2),
        data={'sequence':peptide_string}
    )

linspace = [(i+1)*num_psms for i in range(num_points)]
performance_dict = {}

for tree_type in [TreeType.SORTED_LIST, TreeType.HASHTABLE, TreeType.HASHTABLE_LARGE, TreeType.AVL, TreeType.RB, TreeType.BINARY]:
    performance_dict[tree_type] = {'add_time':[], 'add_time_per_psm':[], 'search_time':[], 'search_time_per_psm':[],
                                   'remove_time':[], 'remove_time_per_psm':[], 'save_time':[], 'save_time_per_psm':[],
                                   'load_time':[], 'load_time_per_psm':[]}
    tree = psm_tree_constructor(tree_type)
    for i in range(num_points):

        psms = [generate_random_psm() for j in range(num_psms)]

        add_start_time = time.time()
        for psm in psms[:n]:
            tree.add(psm)
        add_time = time.time() - add_start_time
        add_time_per_psm = add_time/n

        for psm in psms[n:]:
            tree.add(psm)

        performance_dict[tree_type]['add_time'].append(add_time)
        performance_dict[tree_type]['add_time_per_psm'].append(add_time_per_psm)

        search_start_time = time.time()
        for psm in psms[:n]:
            _ = tree._search(get_mz_bounds(psm.mz, 50),
                            get_rt_bounds(psm.rt, 100),
                            get_ook0_bounds(psm.ook0, 0.05))
        search_time = time.time() - search_start_time
        search_time_per_psm = search_time/n

        performance_dict[tree_type]['search_time'].append(search_time)
        performance_dict[tree_type]['search_time_per_psm'].append(search_time_per_psm)

        remove_start_time = time.time()
        for psm in psms[:n]:
            _ = tree.remove(psm)
        remove_time = time.time() - remove_start_time
        remove_time_per_psm = remove_time/n

        performance_dict[tree_type]['remove_time'].append(remove_time)
        performance_dict[tree_type]['remove_time_per_psm'].append(remove_time_per_psm)

        # add back removed psms
        for psm in psms[:n]:
            tree.add(psm)

        save_start_time = time.time()
        tree.save('tmp.txt')
        save_time = time.time() - save_start_time

        performance_dict[tree_type]['save_time'].append(save_time)
        performance_dict[tree_type]['save_time_per_psm'].append(save_time/n)

        tree.clear()

        load_start_time = time.time()
        tree.load('tmp.txt')
        load_time = time.time() - load_start_time

        performance_dict[tree_type]['load_time'].append(load_time)
        performance_dict[tree_type]['load_time_per_psm'].append(load_time/n)


psms = [generate_random_psm() for j in range(num_psms)]
mz_list = [psm.mz for psm in psms]
ook0_list = [psm.ook0 for psm in psms]

plt.scatter(mz_list,ook0_list)
plt.show()

f, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 1)
f.tight_layout()

for tree_type in performance_dict:
    ax1.plot(linspace, performance_dict[tree_type]['add_time_per_psm'], label=tree_type)

ax1.title.set_text('Add Time per PSM')

for tree_type in performance_dict:
    ax2.plot(linspace, performance_dict[tree_type]['search_time_per_psm'], label=tree_type)
ax2.title.set_text('Search Time per PSM')

for tree_type in performance_dict:
    ax3.plot(linspace, performance_dict[tree_type]['remove_time_per_psm'], label=tree_type)
ax3.title.set_text('Remove Time per PSM')

for tree_type in performance_dict:
    ax4.plot(linspace, performance_dict[tree_type]['save_time'], label=tree_type)
ax4.title.set_text('Save Time')

for tree_type in performance_dict:
    ax5.plot(linspace, performance_dict[tree_type]['load_time'], label=tree_type)
ax5.title.set_text('Load Time')

plt.legend()
plt.show()





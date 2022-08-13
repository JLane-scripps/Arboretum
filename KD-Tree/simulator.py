import os
import time
import random
from arborist import PSMArborist, TreeType, PSM
from boundary import get_mz_bounds, get_rt_bounds, get_ook0_bounds


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


""" -----Variables----- """
PPM = 50
RT_OFF = 250
OOK0_TOL = 5

tree_type = TreeType.SORTED_LIST  # choose which tree type to simulate here
runs = 5  # the number of times we will add [leaf_count] leaves (PSMs) to the tree
leaves = []  # the list of psms (aka 'leaves' on the tree)
leaf_count = 50_000  # The number of PSM's we will generate and add to the tree in each run.
canopy = 0  # Total number of leaves, updated after each run. (leaf_count * runs so far)
performance_dict = {}
out_directory = "leaves"

""" ----- Simulation -----"""
for x in range(runs):
    arb = PSMArborist(tree_type)

    # generate new 50k leaves
    leaves = []
    for i in range(leaf_count):
        leaves.append(generate_random_psm())

    # get add time
    start_time = time.time()
    for psm in leaves:
        arb.add(psm)
    add_time = (time.time() - start_time)
    print(f"N: {canopy:,}, add_time: {add_time:,}")

    # update canopy size
    canopy += leaf_count

    query = [generate_random_psm() for i in range(int(leaf_count * 0.05))]  # 5% of the leaf_count (ex: 50k -> 500)

    # get search time
    start_time = time.time()
    for psm in query:
        query_results = arb.search(psm.charge,
                                   [get_mz_bounds(psm.mz, PPM).lower, get_mz_bounds(psm.mz, PPM).upper],
                                   [get_rt_bounds(psm.rt, RT_OFF).lower, get_rt_bounds(psm.rt, RT_OFF).upper],
                                   [get_ook0_bounds(psm.ook0, OOK0_TOL).lower,
                                    get_ook0_bounds(psm.ook0, OOK0_TOL).upper])

    search_time = (time.time() - start_time)
    print(f"N: {canopy:,}, search_time: {search_time}")
    print(f"N: {canopy:,}, search_time_per_psm: {(time.time() - start_time) / canopy}")

    # save time
    start_time = time.time()
    arb.save(out_directory)
    save_time = (time.time() - start_time)
    print(f"N: {canopy:,}, save time: {save_time:,}")

    # load time
    start_time = time.time()
    arb = PSMArborist(tree_type)
    arb.load(out_directory)
    load_time = (time.time() - start_time)
    print(f"N: {canopy:,}, load time: {load_time:,}")

    performance_dict[canopy] = {'add_time': add_time, 'search_time': search_time, 'save_time': save_time, 'load_time':load_time}

print(performance_dict)

file_name = os.path.join('Sim Log/Sim ' + str(tree_type.name) + '.txt')
with open(file_name, 'w') as newfile:
    newfile.write(str(performance_dict))
newfile.close()

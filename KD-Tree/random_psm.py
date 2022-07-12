import random

def generate_random_psm():
    return {
        'mono_mz': random.uniform(400, 1800),
        'rt': random.uniform(0, 10_000),
        'ook0': random.uniform(0.4, 1.8),
        'charge': random.randint(1, 5),
        'ms2_id': 1
    }
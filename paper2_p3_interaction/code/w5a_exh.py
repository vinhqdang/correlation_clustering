"""Exhaustive check (all set partitions, n=6..9) of the anatomy lemma on OPTIMAL partitions, and (control) on
partitions of cost OPT+1, OPT+2 to see which statements need optimality."""
import sys, random, numpy as np
from w5a_lib import *
def run(n, seeds, fams):
    stats = dict(inst=0, opt_parts=0, viol_opt=0, t1=0, t2=0, tiny=0, ctrl_parts=0, ctrl_viol={}, forcedT2=0, inst_t2=0)
    lab_t, same, P = part_table(n)
    for sd in seeds:
        rng = random.Random(f"w5a-{n}-{sd}")
        A = gen_graph(rng, n, rng.choice(fams))
        if rng.random() < .5: A = noise(A, rng, rng.randint(1, 3))
        H = greedy_extend(A, rand_packing(A, rng, stop=rng.choice([0, .3, .6]), prefer_overlap=rng.random() < .5), rng, prefer_overlap=rng.random() < .5)
        if not H or not is_pair_maximal(A, H): continue
        pk = packed_matrix(n, H)
        lab, same_, cost = all_costs(A); opt = cost.min()
        stats['inst'] += 1
        anyT2free = False; anyopt = False
        for i in np.nonzero(cost <= opt + 2)[0]:
            r = anatomy(A, pk, lab_t[i])
            if cost[i] == opt:
                stats['opt_parts'] += 1; anyopt = True
                stats['t1'] += r['t1']; stats['t2'] += r['t2']; stats['tiny'] += r['tiny_pairs']
                if r['viol']: stats['viol_opt'] += 1; print('VIOL', n, sd, r['viol'][:3], flush=True)
                if r['t2'] == 0: anyT2free = True
            else:
                stats['ctrl_parts'] += 1
                for v in r['viol']: stats['ctrl_viol'][v[0]] = stats['ctrl_viol'].get(v[0], 0) + 1
        if anyopt and not anyT2free: stats['forcedT2'] += 1
    return stats
if __name__ == '__main__':
    n = int(sys.argv[1]); N = int(sys.argv[2]); s0 = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    print(n, run(n, range(s0, s0 + N), ['gnp', 'noisy', 'hub', 'sparse', 'fano']))

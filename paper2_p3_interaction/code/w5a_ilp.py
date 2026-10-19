"""ILP (HiGHS) check of the anatomy lemma on OPTIMAL partitions, n=10..13, structured + random families.
For each instance: several optimal partitions (random side objective, cost <= OPT) -> anatomy().  Also P13 of wave2_search."""
import sys, random, itertools, numpy as np
from w5a_lib import *
from w5a_loc import structured
from a4_lib import ilp, labels_from_z

def optima(A, rng, k=6):
    n = len(A); P = prs(n); r0 = ilp(A); opt = r0[1]; out = [labels_from_z(n, r0[0], r0[2])]
    for _ in range(k):
        obj = rng.uniform(-1, 1, size=len(P)) * 0.01 + (1 - 2 * np.array([A[i][j] for i, j in P], float)) * 1000
        r = ilp(A, obj=obj, cost_le=opt)
        if r: out.append(labels_from_z(n, r[0], r[2]))
    return opt, out

def p13():
    n = 13; A = np.zeros((n, n), np.int8)
    for i in range(n - 1): A[i, i + 1] = A[i + 1, i] = 1
    return A, [(1, 2, 3), (4, 5, 6), (7, 8, 9), (10, 11, 12)]

if __name__ == '__main__':
    mode = sys.argv[1]; n = int(sys.argv[2]); N = int(sys.argv[3])
    inst = parts = viol = t1 = t2 = tiny = 0; cvio = {}
    for sd in range(N):
        rng = random.Random(f"w5ailp-{mode}-{n}-{sd}"); nrng = np.random.default_rng(sd)
        if mode == 'p13':
            A, H = p13(); n = 13
        else:
            A = structured(rng, n, mode) if mode in ('path', 'cycle', 'pathplus', 'caterpillar', 'twopaths') else gen_graph(rng, n, mode)
            if mode not in ('path', 'cycle') and rng.random() < .4: A = noise(A, rng, rng.randint(1, 3))
            H = greedy_extend(A, rand_packing(A, rng, stop=rng.choice([0, .3, .6]), prefer_overlap=rng.random() < .5), rng, prefer_overlap=rng.random() < .5)
        if not H or not is_pair_maximal(A, H): continue
        pk = packed_matrix(n, H); inst += 1
        opt, labs = optima(A, nrng)
        for lab in labs:
            r = anatomy(A, pk, lab); parts += 1; t1 += r['t1']; t2 += r['t2']; tiny += r['tiny_pairs']
            if r['viol']: viol += 1; print('VIOL', mode, n, sd, r['viol'][:3])
        if mode == 'p13': break
    print(mode, n, 'instances', inst, 'optimal partitions', parts, 'violations', viol, 'T1 pairs', t1, 'T2 pairs', t2, 'tiny-tiny cluster pairs', tiny)

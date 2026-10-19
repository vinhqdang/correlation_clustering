"""Smallest counterexamples when 'maximal' is weakened (vertex-disjoint maximal / arbitrary packing).
Exhaustive over all partitions; reports n, edges, H and the optima."""
import random, itertools, numpy as np
import a3a_lib as L
rng = random.Random(7)
best = {}
for n in range(4, 8):
    for trial in range(6000):
        A = L.gen_gnp(rng, n)
        for hyp in ('vertexmax', 'empty_or_partial'):
            if hyp == 'vertexmax': H = L.vertex_greedy_packing(A, rng)
            else:
                H = L.greedy_packing(A, rng); H = H[:rng.randint(0, len(H))]
            if hyp == 'vertexmax' and L.is_pair_maximal(A, H): continue
            VH, R, Q, cls = L.info(A, H)
            if len(R) < 2: continue
            costs = L.costs_all(A); OPT = costs.min(); opt = np.nonzero(costs == OPT)[0]
            S = L.same_matrix(n); pi = L.pidx(n)
            rr = [pi[(r, s)] for r, s in itertools.combinations(R, 2)]
            ed = np.array([A[i, j] for i, j in itertools.combinations(range(n), 2)])
            viol = S[:, rr].astype(int) @ (1 - ed[rr]) + (1 - S[:, rr].astype(int)) @ ed[rr]
            if viol[opt].min() > 0 and hyp not in best:
                E = [(i, j) for i, j in itertools.combinations(range(n), 2) if A[i, j]]
                best[hyp] = (n, E, H, sorted(R), int(OPT), len(H), [tuple(L.rgs(n)[o]) for o in opt][:4])
    print('n', n, {k: v[0] for k, v in best.items()})
    if len(best) == 2: break
for k, v in best.items(): print(k, v)
# P5 example for the c = 1 remark
A = L.sym(5)
for i in range(4): L.add(A, i, i + 1)
H = [(0, 1, 2)]; print('P5: vertex-maximal', not L.is_pair_maximal(A, H), 'pair-extension', L.greedy_packing(A, random.Random(1), None, start=H))

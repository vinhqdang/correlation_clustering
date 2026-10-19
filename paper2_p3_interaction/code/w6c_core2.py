"""Generalised small-core test with ILP (n up to ~16).  Lock = a hard component K (component of the graph of unpacked edges)
that is forced to be a whole cluster (what a twin stub x' adjacent to all of K costs nothing for).  S is a face if some optimal
solution makes all components of S whole clusters.  Report minimal non-faces of size >= 3 (pairwise faces, joint non-face)."""
import sys, itertools, random, numpy as np, networkx as nx
from w6c_lib import *
from w6c_core import max_packing_ilp

def closed_fixed(n, comps):
    fixed = {}
    for K in comps:
        Ks = set(K)
        for a in K:
            for t in range(n):
                if t not in Ks: fixed[key(a, t)] = 0
    return fixed

def test_instance(A, rng, trials=3, maxcomp=None):
    n = len(A)
    opt, _ = ilp_opt(A)
    out = []
    for _ in range(trials):
        H = max_packing_ilp(A, rng)
        if len(H) != opt: return None
        used = set()
        for t in H: used |= pairs_of(t)
        hg = nx.Graph(); hg.add_nodes_from(range(n))
        for i, j in prs(n):
            if A[i][j] and key(i, j) not in used: hg.add_edge(i, j)
        comps = [tuple(sorted(c)) for c in nx.connected_components(hg)]
        # hard components must be cliques for a twin stub to be harmless (else all-hard P3 -> not tight anyway)
        ok1 = []
        for K in comps:
            r = ilp_opt(A, fixed=closed_fixed(n, [K]))
            if r and r[0] == opt: ok1.append(K)
        comp = nx.Graph(); comp.add_nodes_from(ok1)
        for K1, K2 in itertools.combinations(ok1, 2):
            r = ilp_opt(A, fixed=closed_fixed(n, [K1, K2]))
            if r and r[0] == opt: comp.add_edge(K1, K2)
        viol = None
        for cl in nx.find_cliques(comp):
            if len(cl) < 3: continue
            r = ilp_opt(A, fixed=closed_fixed(n, cl))
            if not (r and r[0] == opt): viol = cl; break
        out.append((len(comps), len(ok1), max((len(c) for c in nx.find_cliques(comp)), default=0), viol))
        if viol: return ('VIOLATION', H, comps, viol)
    return out

def planted(rng, n, k, flips):
    sizes = [1] * k
    for _ in range(n - k): sizes[rng.randrange(k)] += 1
    lab = [i for i, s in enumerate(sizes) for _ in range(s)]
    A = np.zeros((n, n), np.int8)
    for i, j in prs(n):
        if lab[i] == lab[j]: A[i, j] = A[j, i] = 1
    for _ in range(flips):
        i, j = rng.sample(range(n), 2); A[i, j] ^= 1; A[j, i] ^= 1
    return A

if __name__ == '__main__':
    n = int(sys.argv[1]); seed = int(sys.argv[2]); T = float(sys.argv[3]); mode = sys.argv[4]
    import time
    rng = random.Random(seed); nrng = np.random.default_rng(seed); t0 = time.time()
    stats = dict(inst=0, tight=0, maxclique=0); 
    while time.time() - t0 < T:
        if mode == 'planted':
            A = planted(rng, n, rng.randint(2, max(2, n // 3)), rng.randint(3, n))
        else:
            g = nx.gnm_random_graph(n, rng.randint(n, int(2.2 * n)), seed=rng.randrange(10**9)); A = np.array(nx.to_numpy_array(g), dtype=np.int8)
        stats['inst'] += 1
        r = test_instance(A, nrng, trials=2)
        if r is None: continue
        stats['tight'] += 1
        if isinstance(r, tuple): print('VIOLATION', r, flush=True); print(A.tolist(), flush=True); break
        stats['maxclique'] = max(stats['maxclique'], max(x[2] for x in r))
    print(n, mode, stats, flush=True)

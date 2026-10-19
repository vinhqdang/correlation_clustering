"""Strengthened HUB: for a YES instance (OPT <= |H|+l), the graph Gamma_u on N_h(u) whose edges are the HARD NON-ADJACENT pairs
has matching number nu <= l (so vertex cover <= 2l), hence alpha(G[N_h(u)]) <= 2l + c + 1.  Also records max tw(W) per w.
usage: python3 a4_hub2.py seed ninst nmin nmax [ilp]"""
import sys, random, collections, itertools
import numpy as np, networkx as nx, a4_lib as L
seed, ninst, nmin, nmax = map(int, sys.argv[1:5]); use_ilp = len(sys.argv) > 5
rng = random.Random(seed); C = collections.Counter(); bad = []; twmax = collections.defaultdict(int); slack = collections.Counter()
for it in range(ninst):
    n = rng.randint(nmin, nmax); A = L.gen_graph(rng, n)
    if rng.random() < .5: A = L.noise(A, rng, rng.randint(1, 3))
    mode = rng.choice(['maximal', 'partial', 'overlap'])
    H = L.rand_packing(A, rng, stop=.4 if mode == 'partial' else 0.0, prefer_overlap=mode == 'overlap')
    if mode != 'partial': H = L.greedy_extend(A, H, rng, mode == 'overlap')
    pk = L.packed_matrix(n, H); h = len(H); c = L.cmax(n, H)
    OPT = L.ilp(A)[1] if use_ilp else int(L.all_costs(A)[2].min())
    if 0 < h <= 13:
        w = L.tw_exact(L.ih_adj(H)); twmax[w] = max(twmax[w], L.tw_exact(L.soft_adj(n, H)))
    for u in range(n):
        Nh = [x for x in range(n) if A[u][x] and not pk[u, x]]
        g = nx.Graph(); g.add_nodes_from(Nh)
        g.add_edges_from((x, y) for x, y in itertools.combinations(Nh, 2) if not A[x][y] and not pk[x, y])
        nu = len(nx.max_weight_matching(g, maxcardinality=True))
        al = L.max_indep(Nh, lambda x, y: A[x][y])
        for l in range(4):
            if OPT <= h + l:
                C['checks'] += 1
                if nu > l: bad.append(('nu', seed, it, u, l, nu))
                if al > 2 * l + c + 1: bad.append(('alpha-additive', seed, it, u, l, al, c))
                slack[(l, 2 * l + c + 1 - al)] += (l == OPT - h)   # only count when l is the exact excess
print(dict(C)); print('max tw(W) by w:', dict(sorted(twmax.items())))
print('VIOLATIONS', len(bad), bad[:5])

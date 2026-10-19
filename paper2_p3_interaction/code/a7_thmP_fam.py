"""Audit 7: expander family of dproute Prop (l=0,w=0,c=1) : weights of cross soft pairs, components of Gamma^l, tw(Gamma')."""
import itertools, random
import numpy as np, networkx as nx
from networkx.algorithms.approximation import treewidth_min_fill_in
import a7_lib as L, a7_thmP as T

def orient_in12(F, rng):
    """orientation with every in-degree in {1,2} (random search + repair by flipping paths)."""
    edges = list(F.edges())
    for _ in range(2000):
        ori = [(u, v) if rng.random() < .5 else (v, u) for u, v in edges]
        for _ in range(500):
            indeg = {v: 0 for v in F}
            for u, v in ori: indeg[v] += 1
            bad = [v for v in F if indeg[v] not in (1, 2)]
            if not bad: return ori
            v = rng.choice(bad)
            inc = [i for i, (a, b) in enumerate(ori) if b == v] if indeg[v] == 3 else [i for i, (a, b) in enumerate(ori) if a == v]
            i = rng.choice(inc); ori[i] = (ori[i][1], ori[i][0])
    raise RuntimeError

def build(F, ori):
    ids = {}; cnt = 0; cl = {v: [] for v in F}; cross = []; H = []
    for k, (p, q) in enumerate(ori):
        a, b, c = cnt, cnt + 1, cnt + 2; cnt += 3
        cl[p].append(a); cl[q] += [b, c]; cross.append((a, b)); H.append((a, b, c))
    A = np.zeros((cnt, cnt), int)
    for v, vs in cl.items():
        for x, y in itertools.combinations(vs, 2): A[x][y] = A[y][x] = 1
    for a, b in cross: A[a][b] = A[b][a] = 1
    return A, H, cl

for name, F in (("K4", nx.complete_graph(4)), ("Petersen", nx.petersen_graph()), ("cubic20", nx.random_regular_graph(3, 20, seed=3)), ("cubic40", nx.random_regular_graph(3, 40, seed=4))):
    rng = random.Random(1); ori = orient_in12(F, rng)
    A, H, cl = build(F, ori)
    n = len(A)
    assert L.is_packing(A, H) and L.is_pair_maximal(A, H)
    S = L.soft_of(n, H); Dl, Gm = L.delta_gamma(A, S)
    inK = {x: v for v, vs in cl.items() for x in vs}
    w = sorted(set(len(Dl[p]) for p in L.prs(n) if S[p] and inK[p[0]] != inK[p[1]]))
    for l in (0, 1, 2, 3, 4, 5):
        st = T.setup(A, H, l)
        comps = [len(c) for c in nx.connected_components(st['GL'])]
        print(name, "n", n, "D", st['D'], "BIG", len(st['BIG']), "cross soft weights", w, "l", l, "max comp of Gamma^l", max(comps), "tw(Gamma^l) <=", treewidth_min_fill_in(st['GL'])[0])
    # Gamma' treewidth upper bound & quotient contains F
    G2 = nx.Graph(); G2.add_nodes_from(range(n))
    for i, j in L.prs(n):
        if A[i][j] or S[i, j]: G2.add_edge(i, j)
    Q = nx.quotient_graph(G2, [set(vs) for vs in cl.values()])
    print("  quotient edges == F edges:", Q.number_of_edges() == F.number_of_edges(), " tw(Gamma') <= (min-fill)", treewidth_min_fill_in(G2)[0], " tw(F) <= ", treewidth_min_fill_in(F)[0])
    print("  OPT==|H| check via partition cost:", L.opt_subset(A) if n <= 12 else "skip", len(H))

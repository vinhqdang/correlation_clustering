"""Spine family Sp(F) (paper Prop. five / app_spine.tex), built from scratch; checks (S1) and the one-exchange repack
H' = H - X_0 + Y + Z  with Y=(y,u_0,b_0), Z=(b_0,u_1,c_1):  |H'| = |H|+1 = OPT  =>  excess 0, tw(I(H')) small."""
import random, itertools, sys, time
import networkx as nx
from w8a_lib import *
def orient(F, rng):
    # in-degree <= 2 orientation via bipartite matching/flow
    D = nx.DiGraph()
    for i, (u, v) in enumerate(F.edges()):
        D.add_edge('s', ('e', i), capacity=1)
        D.add_edge(('e', i), ('v', u), capacity=1); D.add_edge(('e', i), ('v', v), capacity=1)
    for v in F.nodes(): D.add_edge(('v', v), 't', capacity=2)
    val, fl = nx.maximum_flow(D, 's', 't')
    assert val == F.number_of_edges()
    out = []
    for i, (u, v) in enumerate(F.edges()):
        head = u if fl[('e', i)][('v', u)] > 0 else v
        tail = v if head == u else u
        out.append((tail, head))
    return out
def build(F, rng):
    arcs = orient(F, rng); t0 = arcs[0][0]
    for _ in range(2000):
        rest = arcs[:]; rng.shuffle(rest)
        arcs_all = [(t0, '*')] + rest
        if all(arcs_all[i][1] != arcs_all[i + 1][1] for i in range(len(arcs_all) - 1)) and arcs_all[1][1] != t0: break
    else: return None
    N = len(arcs_all) - 1
    idx = {}; nv = 0
    for i in range(N + 1):
        for nm in 'abcu': idx[(nm, i)] = nv; nv += 1
    z = nv; y = nv + 1; n = nv + 2
    K = {}
    for i, (t, h) in enumerate(arcs_all):
        K.setdefault(t, []).append(idx[('a', i)])
        for nm in 'bcu': K.setdefault(h, []).append(idx[(nm, i)])
    E = set()
    for v, vs in K.items():
        for x, w in itertools.combinations(vs, 2): E.add(key(x, w))
    H = []
    for i in range(N + 1):
        a, b, c, u = (idx[(nm, i)] for nm in 'abcu')
        E.add(key(a, b))
        un = idx[('u', i + 1)] if i < N else z
        E.add(key(b, un))
        H.append((a, b, c)); H.append((u, b, un))
    E.add(key(y, idx[('u', 0)])); E.add(key(y, idx[('c', 0)]))
    return n, E, H, idx, y, K, arcs_all
def run(name, F, rng):
    r = build(F, rng)
    if r is None: print(name, 'ordering failed'); return
    n, E, H, idx, y, K, arcs = r
    g = G8(n, list(E))
    ok = valid_packing(g, H); pm = pair_maximal(g, H)
    I = Igraph(H)
    # explicit partition of cost |H|+1: cliques K_v+ (v in F'), z alone, y joins K_*+
    lab = {}
    for j, (v, vs) in enumerate(K.items()):
        for x in vs: lab[x] = j
    lab[n - 2] = max(lab.values()) + 1; lab[y] = lab[idx[('u', 0)]]
    part = {}
    for v, l in lab.items(): part.setdefault(l, []).append(v)
    cost = len(edits_of_partition(g, list(part.values())))
    # repack
    u0, b0, c0 = idx[('u', 0)], idx[('b', 0)], idx[('c', 0)]; u1, c1 = idx[('u', 1)], idx[('c', 1)]
    X0 = (u0, b0, idx[('u', 1)]); assert X0 in H
    Y = (y, u0, b0); Z = (b0, u1, c1)
    H2 = [t for t in H if t != X0] + [Y, Z]
    ok2 = valid_packing(g, H2); pm2 = pair_maximal(g, H2)
    I2 = Igraph(H2)
    from networkx.algorithms.approximation import treewidth_min_fill_in
    tw1 = treewidth_min_fill_in(I)[0]; tw2 = treewidth_min_fill_in(I2)[0]
    print('%-14s n=%3d |H|=%3d valid=%s pairmax=%s tw(I(H))<=%d | explicit cost=%d (=|H|+1: %s) | H2 valid=%s |H2|=%d (=|H|+1: %s) pairmax(H2)=%s tw(I(H2))<=%d  => OPT=|H2|=%d, excess 0'
          % (name, n, len(H), ok, pm, tw1, cost, cost == len(H) + 1, ok2, len(H2), len(H2) == len(H) + 1, pm2, tw2, len(H2)))
if __name__ == '__main__':
    rng = random.Random(7)
    run('triangle', nx.complete_graph(3), rng)
    run('K4', nx.complete_graph(4), rng)
    run('K33', nx.complete_bipartite_graph(3, 3), rng)
    run('Petersen', nx.petersen_graph(), rng)
    for nn in (10, 16, 24, 40):
        for s in range(2):
            F = nx.random_regular_graph(3, nn, seed=s)
            run('cubic%d_%d' % (nn, s), F, rng)

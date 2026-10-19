"""Audit 6, Prop 9 (far-default-link family). Own construction from the verbal description.
Part 1: exhaustive (all set partitions) on small F (n = 3|E(F)| <= 10 -> |E(F)| <= 3): OPT = |H|, H induced, pair-disjoint, pair-maximal.
Part 2: structural checks on cubic F (K4, K33, Petersen, random cubic): induced, pair-disjoint, pair-maximal, cost({K_v}) = |H|,
        packing lower bound => OPT = |H|; max |K_v|, max hard degree, B = 2D+3l ; F is a minor of Gamma' (explicit contraction check) ;
        and whether the family trivialises under unit propagation (every K_v is a hard component, no edited hard pair possible).
usage: python3 -I a6_p9.py"""
import sys, itertools, random
import numpy as np, networkx as nx
sys.path.insert(0, '.')
from a6_lib import *

def build(F, orient=None):
    verts = {}; cnt = 0; trip = []; K = {v: [] for v in F}
    for (u, v) in F.edges():
        if orient and orient.get((u, v)) == 'rev': u, v = v, u
        a, b, c = cnt, cnt + 1, cnt + 2; cnt += 3; trip.append((a, b, c)); K[u].append(a); K[v] += [b, c]
    A = np.zeros((cnt, cnt), int)
    for v, mem in K.items():
        for x, y in itertools.combinations(mem, 2): A[x][y] = A[y][x] = 1
    for a, b, c in trip: A[a][b] = A[b][a] = 1
    return A, trip, K

def checks(A, trip, K):
    n = len(A); S = soft_matrix(n, trip)
    ind = all(A[a][b] and A[b][c] and not A[a][c] for a, b, c in trip)
    disj = len({p for t in trip for p in p3pairs(t)}) == 3 * len(trip)
    pm = all(any(S[x, y] for x, y in ((a, b), (b, c), (a, c))) for a, b, c in induced_p3s(A))
    lab = {x: v for v, mem in K.items() for x in mem}
    cost = sum(1 for x, y in pairs(n) if (lab[x] == lab[y]) != bool(A[x][y]))
    return S, ind, disj, pm, cost

if __name__ == '__main__':
    ok = True
    for name, F in (('triangle', nx.cycle_graph(3)), ('star3', nx.star_graph(3)), ('path4', nx.path_graph(4)), ('2+1', nx.Graph([(0, 1), (2, 3), (3, 4)]))):
        if 3 * F.number_of_edges() > 10: continue
        for rev in (False, True):
            A, trip, K = build(F, {e: 'rev' for e in F.edges()} if rev else None)
            S, ind, disj, pm, cost = checks(A, trip, K)
            P = partitions(len(A)); opt = int(part_costs(A, P).min())
            print(name, 'n=%d' % len(A), 'induced', ind, 'disjoint', disj, 'pair-max', pm, 'cost(K_v)', cost, '|H|', len(trip), 'OPT(brute)', opt)
            ok &= ind and disj and pm and cost == len(trip) == opt
    for name, F in (('K4', nx.complete_graph(4)), ('K33', nx.complete_bipartite_graph(3, 3)), ('Petersen', nx.petersen_graph()),
                    ('cubic20', nx.random_regular_graph(3, 20, seed=2)), ('cubic40', nx.random_regular_graph(3, 40, seed=4))):
        A, trip, K = build(F)
        S, ind, disj, pm, cost = checks(A, trip, K); n = len(A)
        comp, hg = hard_components(A, S); dh = hg.sum(1); D = int(S.sum(1).max())
        # every K_v is exactly one hard component and clean
        comps = {}
        for v in range(n): comps.setdefault(comp[v], set()).add(v)
        same = sorted(map(sorted, comps.values())) == sorted(map(sorted, K.values()))
        clean = all(A[x][y] or S[x, y] for C in comps.values() for x, y in itertools.combinations(sorted(C), 2))
        # Gamma' and the minor claim
        G = nx.Graph(); G.add_nodes_from(range(n))
        for x, y in pairs(n):
            if A[x][y] or S[x, y]: G.add_edge(x, y)
        Q = nx.quotient_graph(G, [set(m) for m in K.values()], relabel=False)
        qe = {frozenset(map(frozenset, e)) for e in Q.edges()}
        Fe = {frozenset((frozenset(K[u]), frozenset(K[v]))) for u, v in F.edges()}
        print(name, 'n=%d' % n, 'induced', ind, 'disj', disj, 'pair-max', pm, 'cost', cost, '|H|', len(trip), 'D', D, 'max|K_v|', max(len(m) for m in K.values()),
              'max d_h', int(dh.max()), 'B(l=0)=2D=', 2 * D, 'K_v == hard comps', same, 'clean', clean, 'quotient(Gamma\') == F', qe == Fe)
        ok &= ind and disj and pm and cost == len(trip) and same and clean and qe == Fe
    print('ALL OK' if ok else 'SOMETHING FAILED')

"""Audit 6, Thm 8: (a) literal check of the proof's claim 'L_B subset V'' (b) exceptional clusters vs BIG, (c) all non-BIG clusters of every
optimal partition are Gamma'-connected and <= Z, BIG components are whole clusters, (d) the formula  |E|+sum f(K)+m <= |H|+l  <=> OPT<=|H|+l,
with m computed (i) by brute force over restricted partitions and (ii) by an actual tree-decomposition DP (introduce / forget / JOIN),
(e) DP without the join cross term x1*x2 (what the note literally states) must FAIL somewhere.
usage: python3 -I a6_t8.py N SEED NMIN NMAX"""
import sys, random, itertools
import numpy as np, networkx as nx
from networkx.algorithms.approximation import treewidth_min_fill_in
sys.path.insert(0, '.')
from a6_lib import *

def fC(A, C):
    k = len(C); e = sum(int(A[x][y]) for x, y in itertools.combinations(C, 2)); return k * (k - 1) // 2 - 2 * e

def connected(G, C):
    C = set(C); s = next(iter(C)); seen = {s}; st = [s]
    while st:
        v = st.pop()
        for u in G[v]:
            if u in C and u not in seen: seen.add(u); st.append(u)
    return seen == C

def brute_m(A, G, V, Z, connected_only):
    """min sum f(C) over partitions of V with clusters <= Z (and Gamma'-connected if requested)"""
    V = list(V); best = [10 ** 9]
    if not V: return 0
    def rec(i, blocks, cost_so_far):
        if i == len(V):
            tot = 0
            for C in blocks:
                if len(C) > Z: return
                if connected_only and not connected(G, C): return
                tot += fC(A, C)
            best[0] = min(best[0], tot); return
        v = V[i]
        for b in blocks:
            if len(b) < Z: b.append(v); rec(i + 1, blocks, 0); b.pop()
        blocks.append([v]); rec(i + 1, blocks, 0); blocks.pop()
    rec(0, [], 0); return best[0]

def dp_m(A, G, Z, join_term=True):
    """tree-decomposition DP for min sum f(C), clusters <= Z; per connected component of G. States: (block label tuple over sorted bag, sizes)."""
    total = 0
    for comp in nx.connected_components(G):
        g = G.subgraph(comp).copy()
        if g.number_of_nodes() == 1: continue   # singleton: f = 0
        w, T = treewidth_min_fill_in(g)
        root = next(iter(T.nodes))
        def canon(bag, lab, size):   # lab: dict v->block id ; size: dict block->size ; canonical relabel
            order = sorted(bag); mp = {}; L = []
            for v in order:
                b = lab[v]
                if b not in mp: mp[b] = len(mp)
                L.append(mp[b])
            sz = [0] * len(mp)
            for b, i in mp.items(): sz[i] = size[b]
            return (tuple(L), tuple(sz))
        def decode(bag, st):
            order = sorted(bag); L, sz = st
            return {v: L[i] for i, v in enumerate(order)}, {i: s for i, s in enumerate(sz)}
        def forget(table, bag_old, bag_new):
            out = {}
            for st, cst in table.items():
                lab, size = decode(bag_old, st)
                keep = {v: lab[v] for v in bag_new}
                sz = {b: size[b] for b in set(keep.values())}
                k = canon(bag_new, keep, sz)
                if k not in out or out[k] > cst: out[k] = cst
            return out
        def introduce(table, bag, v):
            out = {}
            nb = sorted(bag); bag2 = set(bag) | {v}
            for st, cst in table.items():
                lab, size = decode(bag, st)
                # new block
                l2 = dict(lab); l2[v] = 99; s2 = dict(size); s2[99] = 1
                k = canon(bag2, l2, s2)
                if k not in out or out[k] > cst: out[k] = cst
                for b in set(lab.values()):
                    s = size[b]
                    if s + 1 > Z: continue
                    kk = sum(1 for u in bag if lab[u] == b and A[u][v])
                    l2 = dict(lab); l2[v] = b; s2 = dict(size); s2[b] = s + 1
                    k = canon(bag2, l2, s2); c2 = cst + s - 2 * kk
                    if k not in out or out[k] > c2: out[k] = c2
            return out
        def solve(t, parent):
            bag = set(t)
            tab = {((), ()): 0}; cur = set()
            kids = [c for c in T.neighbors(t) if c != parent]
            tabs = []
            for c in kids:
                ct = solve(c, t); cbag = set(c)
                ct = forget(ct, cbag, cbag & bag); cur_b = cbag & bag
                for v in sorted(bag - cbag): ct = introduce(ct, cur_b, v); cur_b = cur_b | {v}
                tabs.append(ct)
            if not tabs:
                cur_b = set(); ct = {((), ()): 0}
                for v in sorted(bag): ct = introduce(ct, cur_b, v); cur_b = cur_b | {v}
                return ct
            res = tabs[0]
            for ct in tabs[1:]:
                new = {}
                for s1, c1 in res.items():
                    for s2, c2 in ct.items():
                        if s1[0] != s2[0]: continue
                        L = s1[0]; nbk = len(s1[1]); ok = True; sz = []; cost = c1 + c2
                        order = sorted(bag)
                        for i in range(nbk):
                            mem = [order[j] for j in range(len(order)) if L[j] == i]; b = len(mem)
                            x1, x2 = s1[1][i] - b, s2[1][i] - b
                            ssum = b + x1 + x2
                            if ssum > Z: ok = False; break
                            sz.append(ssum)
                            fb = fC(A, mem)
                            cost -= fb
                            if join_term: cost += x1 * x2
                        if not ok: continue
                        k = (L, tuple(sz))
                        if k not in new or new[k] > cost: new[k] = cost
                res = new
            return res
        tab = solve(root, None)
        total += min(tab.values())
    return total

def run(ninst, seed, nmin, nmax, dpmax=9):
    rng = random.Random(seed)
    st = dict(inst=0, claim_LB_in_V_false=0, exc_meets_big=0, big_nontriv=0, eq_checks=0, dp_checks=0, dp_nojoin_wrong=0, fail=[], maxslack_Z=-99)
    for it in range(ninst):
        n = rng.randint(nmin, nmax); A = rand_graph(n, rng)
        if A.sum() == 0: continue
        H = greedy_pair_maximal(A, rng)
        if not H: continue
        S = soft_matrix(n, H); D = int(S.sum(1).max()); comp, hg = hard_components(A, S); dh = hg.sum(1)
        Pt = partitions(n); c = part_costs(A, Pt); opt = int(c.min()); hp, masks = hard_edit_masks(A, S, Pt)
        nE = int(A.sum() // 2); st['inst'] += 1
        for l in range(max(0, opt - len(H)), opt - len(H) + 2):
            B = 2 * D + 3 * l; Z = 3 * D + 3 * l + 1
            bigc = {comp[v] for v in range(n) if dh[v] > B}
            BIG = [v for v in range(n) if comp[v] in bigc]; Vp = [v for v in range(n) if comp[v] not in bigc]
            LB = [v for v in range(n) if dh[v] <= B]
            if set(LB) - set(Vp): st['claim_LB_in_V_false'] += 1
            G = nx.Graph(); G.add_nodes_from(Vp)
            for x, y in itertools.combinations(Vp, 2):
                if A[x][y] or S[x, y]: G.add_edge(x, y)
            clean = all(not (not A[x][y] and not S[x, y]) for k in bigc for x, y in itertools.combinations([v for v in range(n) if comp[v] == k], 2))
            sumfK = sum(fC(A, [v for v in range(n) if comp[v] == k]) for k in bigc)
            # (b),(c) on every optimal partition (only meaningful when YES for this l)
            if opt <= len(H) + l:
                if not clean: st['fail'].append(('notclean_but_YES', it))
                for ii in np.nonzero(c == opt)[0]:
                    lab = Pt[ii]; cl = {}
                    for v in range(n): cl.setdefault(int(lab[v]), []).append(v)
                    U = [hp[k] for k in range(len(hp)) if (int(masks[ii]) >> k) & 1]
                    ends = {v for p in U for v in p}
                    for C in cl.values():
                        Cs = set(C)
                        if Cs & set(BIG):
                            st['big_nontriv'] += 1
                            K = {v for v in range(n) if comp[v] in bigc and comp[v] == comp[next(iter(Cs & set(BIG)))]}
                            if Cs != K: st['fail'].append(('BIGcluster', it))
                            continue
                        if Cs & ends and Cs & set(BIG): st['exc_meets_big'] += 1
                        if len(C) > 1 and not connected(G, C): st['fail'].append(('Gconn', it, C))
                        st['maxslack_Z'] = max(st['maxslack_Z'], len(C) - Z)
                        if len(C) > Z: st['fail'].append(('Z', it))
            # (d)
            if clean:
                m1 = brute_m(A, G, Vp, Z, True) if len(Vp) <= 10 else None
                lhs = nE + sumfK + m1 <= len(H) + l
                if lhs != (opt <= len(H) + l): st['fail'].append(('EQ', it, l, lhs, opt))
                st['eq_checks'] += 1
                if n <= dpmax and len(Vp) > 0:
                    mdp = dp_m(A, G, Z, True); mall = brute_m(A, G, Vp, Z, False)
                    if not (mall <= mdp <= m1): st['fail'].append(('DP', it, mall, mdp, m1))
                    st['dp_checks'] += 1
                    mnj = dp_m(A, G, Z, False)
                    if mnj != mdp: st['dp_nojoin_wrong'] += 1
            else:
                if opt <= len(H) + l: st['fail'].append(('x', it))
    return st

if __name__ == '__main__':
    N, seed, a, b = map(int, sys.argv[1:5]); st = run(N, seed, a, b); f = st.pop('fail'); print(st); print('FAILS', len(f), f[:8])

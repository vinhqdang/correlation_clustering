"""w6b_prune.py -- the pruning rule behind Theorem P (notes/wave6_pointers.md).
 (F1)  a soft pair (u,v) with |Delta(u,v)| > l is never inside one cluster of a solution with |S| <= |H| + l.
 Part 1: brute force over all partitions with cost <= |H|+l (l = 0..3) of random pair-maximal instances, n <= 8.
 Part 2: the cubic-graph family of Proposition 'expander' (dproute.tex): for K_4 and the Petersen graph, compute
         |Delta| of every cross soft pair, prune with l = 3, and report the largest connected component of
         Gamma^l = (hard edges + unpruned soft pairs): its size bounds tw(Gamma^l)+1, while tw(Gamma') >= tw(F).
Usage: python3 w6b_prune.py [n] [n_inst] [seed]
"""
import sys, random, itertools
import numpy as np
from w6b_lib import *

def part1(n, ninst, seed):
    rng = random.Random(seed); viol = 0; chk = 0; pruned_pairs = 0; tot_pairs = 0
    for _ in range(ninst):
        adj = rand_graph(n, rng)
        H = random_pm_packing(adj, rng, maxc=rng.choice([1, 2, 3, None]), drop=rng.choice([0, 0.3]))
        if not H: continue
        inst = Inst(adj, H); inst.solve()
        soft = sorted(tuple(sorted(p)) for p in inst.soft)
        w = {pr: len(delta(inst, *pr)) for pr in soft}
        for l in range(0, 4):
            sel = inst.P[inst.cost <= len(H) + l]
            forb = [pr for pr in soft if w[pr] > l]
            tot_pairs += len(soft); pruned_pairs += len(forb)
            for lab in sel:
                chk += 1
                for (u, v) in forb:
                    if lab[u] == lab[v]: viol += 1
    return chk, viol, pruned_pairs, tot_pairs

def expander_family(Fedges, nv):
    """Proposition expander: orient edges of F; for e=(p->q): a_e in K_p, b_e,c_e in K_q; cross edge a_e b_e."""
    verts = {}; cnt = 0
    def newv(name):
        nonlocal cnt
        verts[name] = cnt; cnt += 1
    for i, (p, q) in enumerate(Fedges):
        for t in 'abc': newv((t, i))
    K = {v: [] for v in range(nv)}
    for i, (p, q) in enumerate(Fedges):
        K[p].append(verts[('a', i)]); K[q] += [verts[('b', i)], verts[('c', i)]]
    E = set()
    for v in range(nv):
        for x, y in itertools.combinations(K[v], 2): E.add(frozenset((x, y)))
    for i in range(len(Fedges)): E.add(frozenset((verts[('a', i)], verts[('b', i)])))
    H = [(verts[('a', i)], verts[('b', i)], verts[('c', i)]) for i in range(len(Fedges))]
    return cnt, E, H, K

def part2():
    out = []
    graphs = {'K4': (4, [(0, 1), (1, 2), (2, 3), (3, 0), (0, 2), (1, 3)]),
              'Petersen': (10, [(i, (i + 1) % 5) for i in range(5)] + [(i, i + 5) for i in range(5)] + [(5 + i, 5 + (i + 2) % 5) for i in range(5)])}
    for name, (nv, Fe) in graphs.items():
        n, E, H, K = expander_family(Fe, nv)
        A = np.zeros((n, n), dtype=np.int8)
        for e in E:
            a, b = tuple(e); A[a, b] = A[b, a] = 1
        assert is_packing(A, H) and pair_maximal(A, H)
        inst = Inst(A, H)
        soft = sorted(tuple(sorted(p)) for p in inst.soft)
        wts = {pr: len(delta(inst, *pr)) for pr in soft}
        cross = [pr for pr in soft if not any(set(pr) <= set(K[v]) for v in K)]
        mn = min(wts[pr] for pr in cross)
        l = 3
        keep = [pr for pr in soft if wts[pr] <= l]
        # Gamma^l: hard edges + unpruned soft pairs
        par = list(range(n))
        def f(x):
            while par[x] != x: par[x] = par[par[x]]; x = par[x]
            return x
        for e in E:
            a, b = tuple(e)
            if frozenset((a, b)) not in inst.soft: par[f(a)] = f(b)
        for (a, b) in keep: par[f(a)] = f(b)
        comp = {}
        for v in range(n): comp[f(v)] = comp.get(f(v), 0) + 1
        out.append((name, n, len(H), 'min weight of a cross soft pair', mn, 'l=%d: largest component of Gamma^l' % l, max(comp.values())))
    return out

def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    ninst = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    chk, viol, pp, tp = part1(n, ninst, seed)
    print('part1 n=%d: partitions checked %d, merged pruned pairs: %d violations; pruned %d of %d (pair,l) cases' % (n, chk, viol, pp, tp))
    for r in part2(): print('part2', r)

if __name__ == '__main__':
    main()

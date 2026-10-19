"""w6b_merge.py -- 'component-link consistency'.  In an optimum, a cluster C containing two pieces A != B (pieces = components of kept hard
edges) forces every soft pair between A and B to be MERGED (it lies inside C).  If two such cross soft pairs lie in members at
I(H)-distance >= K (infinite counts), their states must be EQUAL although no bag contains both members: a far equality constraint
between two soft-pair states (the two-sided T2 / insertion-only phenomenon).  Search: instances (n <= 10) in which EVERY optimum has a
cluster with two pieces joined by two cross soft pairs whose members are at I(H)-distance >= K.
Usage: python3 w6b_merge.py n n_inst seed K
"""
import sys, random
from collections import Counter
import numpy as np
from w6b_lib import *

def pieces_of(inst, lab):
    n = inst.n; par = list(range(n))
    def f(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    for a, b in pair_list(n):
        if inst.adj[a][b] and not inst.is_soft(a, b) and lab[a] == lab[b]: par[f(a)] = f(b)
    return [f(v) for v in range(n)]

def member_of_pair(inst, a, b):
    for i, t in enumerate(inst.H):
        if a in t and b in t: return i

def far_merge(inst, lab, K):
    pc = pieces_of(inst, lab); n = inst.n
    if not hasattr(inst, '_md'): inst._md = inst.member_dist()
    res = []
    for a, b in pair_list(n):
        pass
    cross = {}
    for a, b in pair_list(n):
        if lab[a] == lab[b] and pc[a] != pc[b] and inst.is_soft(a, b):
            key = tuple(sorted((pc[a], pc[b])))
            cross.setdefault(key, []).append(member_of_pair(inst, a, b))
    for key, ms in cross.items():
        for i in range(len(ms)):
            for j in range(i + 1, len(ms)):
                if inst._md[ms[i]][ms[j]] >= K: res.append(key)
    return res

def main():
    n, ninst, seed, K = [int(x) for x in sys.argv[1:5]]
    rng = random.Random(seed); st = Counter(); ex = []
    for _ in range(ninst):
        adj = rand_graph(n, rng)
        H = random_pm_packing(adj, rng, maxc=rng.choice([1, 2, 3, None]), drop=rng.choice([0, 0.3]))
        if not H or not pair_maximal(adj, H): continue
        inst = Inst(adj, H); inst.solve(); opts = inst.opt_partitions(); st['inst'] += 1
        flags = [bool(far_merge(inst, lab, K)) for lab in opts]
        if any(flags): st['some optimum has far-linked merge'] += 1
        if all(flags):
            st['EVERY optimum has far-linked merge'] += 1; ex.append((adj, H))
    print('n=%d K=%d' % (n, K), dict(st))
    for adj, H in ex[:3]:
        print('  edges', [(i, j) for i in range(n) for j in range(i + 1, n) if adj[i][j]], 'H', H)

if __name__ == '__main__':
    main()

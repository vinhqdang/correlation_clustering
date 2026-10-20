"""Graph-level search: ring of N units; unit = K4 {c0..c3} (all 6 edges) + slot vertices w_{i-1} (adjacent to subset Sb of the K4) and w_i (adjacent to subset Sa),
+ anchor vertex x_i adjacent to a subset Sx of {w_{i-1}, w_i, c0..c3}.  Packing H = a maximum pair-disjoint induced-P3 packing (ILP, several random ones);
tight if OPT = |H|.  Anchors = x_i if all its edges are packed.  Report rings with a lock core of size N (full anchor set not a face, every N-1 subset a face)."""
import sys, itertools, random, time, numpy as np
from w6c_lib import *
from w6c_core import max_packing_ilp
from w6c_ringK4 import face
def build(N, Sa, Sb, Sx):
    idx = {}
    def V(name):
        if name not in idx: idx[name] = len(idx)
        return idx[name]
    edges = set()
    for i in range(N):
        c = [V(('c', i, t)) for t in range(4)]
        for a, b in itertools.combinations(c, 2): edges.add(key(a, b))
        wl, wr, x = V(('w', i)), V(('w', (i + 1) % N)), V(('x', i))
        for t in Sb: edges.add(key(wl, c[t]))
        for t in Sa: edges.add(key(wr, c[t]))
        for s in Sx:
            tgt = {'wl': wl, 'wr': wr}.get(s) if isinstance(s, str) else c[s]
            edges.add(key(x, tgt))
    A = adj_from_edges(len(idx), edges)
    return A, [idx[('x', i)] for i in range(N)]
def subsets(k, sizes):
    for r in sizes:
        yield from itertools.combinations(range(k), r)
if __name__ == '__main__':
    N = int(sys.argv[1]); part = int(sys.argv[2]); nparts = int(sys.argv[3]); T = float(sys.argv[4])
    rng = np.random.default_rng(part); t0 = time.time(); cnt = 0; tight = 0
    cfgs = [(Sa, Sb, Sx) for Sa in subsets(4, (1, 2, 3)) for Sb in subsets(4, (1, 2, 3)) for Sx in
            ([('wl',), ('wr',), ('wl', 'wr')] + [(c,) for c in range(4)] + [('wl', c) for c in range(4)] + [('wr', c) for c in range(4)])]
    rng.shuffle(cfgs)
    for ci, (Sa, Sb, Sx) in enumerate(cfgs):
        if ci % nparts != part: continue
        if time.time() - t0 > T: break
        A, xs = build(N, Sa, Sb, Sx); cnt += 1
        o = ilp_opt(A)[0]
        H = max_packing_ilp(A, rng)
        if len(H) != o: continue
        tight += 1
        for trial in range(3):
            H = max_packing_ilp(A, rng) if trial else H
            used = set()
            for t in H: used |= pairs_of(t)
            n = len(A)
            anch = [x for x in xs if all(key(x, u) in used for u in range(n) if A[x][u])]
            if len(anch) != N: continue
            def fc(S):
                fixed = {}
                for a in S:
                    for t in range(n):
                        if t != a: fixed[key(a, t)] = 0
                r = ilp_opt(A, fixed=fixed); return r is not None and r[0] == o
            sub = [fc([xs[i] for i in range(N) if i != j]) for j in range(N)]
            if not any(sub): break
            full = fc(xs)
            print('cfg', (Sa, Sb, Sx), 'tight; (N-1)-faces', sub, 'full face', full, flush=True)
            if all(sub) and not full: print('CORE FOUND', (Sa, Sb, Sx), flush=True)
            break
    print('part', part, 'graphs', cnt, 'tight', tight, flush=True)

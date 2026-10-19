"""Ring of N copies of the 5-member flip gadget G* with a 1-member wire between gadgets (new vertex t_j).
G*: K4 {0,1,4,5} (hard edges 04,15), a=3, b=6, anchor x=7; members (6,3,7),(4,1,3),(5,0,3),(4,5,6),(6,0,1).
Wire j joins a_j and b_{j+1}: member on {a_j, t_j, b_{j+1}} with t_j at position `tp`, edges chosen so it is an induced P3 in G."""
import sys, itertools, numpy as np
from w6c_lib import *
from w6c_ringK4 import face
def build(N, wire_shape, extra_edges=()):
    idx = {}
    def v(name):
        if name not in idx: idx[name] = len(idx)
        return idx[name]
    members = []; edges = set()
    for j in range(N):
        g = lambda t: v((t, j))
        K0, K1, K4, K5, a, b, x = (g(t) for t in ('0', '1', '4', '5', 'a', 'b', 'x'))
        for p in ((K0, K1), (K0, K4), (K0, K5), (K1, K4), (K1, K5), (K4, K5), (a, K0), (a, K1), (b, K0), (b, K5), (b, a), (a, x)): edges.add(key(*p))
        members += [(b, a, x), (K4, K1, a), (K5, K0, a), (K4, K5, b), (b, K0, K1)]
    for j in range(N):
        a = v(('a', j)); b = v(('b', (j + 1) % N)); t = v(('t', j))
        W = {'ctr': (a, t, b), 'a_end': (t, a, b), 'b_end': (a, b, t)}[wire_shape]
        members.append(W); edges |= {key(W[0], W[1]), key(W[1], W[2])}
    n = len(idx); A = adj_from_edges(n, edges)
    return A, members, idx
if __name__ == '__main__':
    import itertools
    for shape in ('ctr', 'a_end', 'b_end'):
        for N in (3, 4, 5):
            A, H, idx = build(N, shape)
            ok, why = check_packing(A, H)
            if not ok: print(shape, N, 'invalid', why); continue
            opt = ilp_opt(A)[0]; xs = [idx[('x', j)] for j in range(N)]
            fam = [S for r in range(N + 1) for S in itertools.combinations(range(N), r) if face(A, [xs[i] for i in S], opt)]
            maxf = [S for S in fam if not any(set(S) < set(T) for T in fam)]
            print(shape, 'N', N, 'n', len(A), '|H|', len(H), 'OPT', opt, 'maximal faces', maxf, flush=True)

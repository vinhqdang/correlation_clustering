"""Odd ring of K2 'pieces' Q_i={a_i,b_i} (consecutive pieces form K4's) with a satellite anchor x_i per piece.
Members: S_i=(x_i,a_i,b_i)  [x_i a_i, a_i b_i edges],  M_a^i=(b_{i-1},a_i,a_{i+1}),  M_b^i=(a_{i-1},b_i,b_{i+1}).
Idea: all anchors locked => pieces must be perfectly matched around the ring => odd N infeasible; releasing any one anchor lets Q_i merge with x_i."""
import sys, itertools, numpy as np
from w6c_lib import *
from w6c_ringK4 import face
def build(N):
    idx = {}
    def v(name):
        if name not in idx: idx[name] = len(idx)
        return idx[name]
    a = lambda i: v(('a', i % N)); b = lambda i: v(('b', i % N)); x = lambda i: v(('x', i % N))
    members = []; edges = set()
    for i in range(N):
        members += [(x(i), a(i), b(i)), (b(i - 1), a(i), a(i + 1)), (a(i - 1), b(i), b(i + 1))]
        edges |= {key(x(i), a(i)), key(a(i), b(i))}
        for p in (a(i), b(i)):
            for q in (a(i + 1), b(i + 1)): edges.add(key(p, q))
    n = len(idx); A = adj_from_edges(n, edges)
    return A, members, [idx[('x', i)] for i in range(N)]
if __name__ == '__main__':
    for N in map(int, sys.argv[1:]):
        A, H, xs = build(N)
        ok, why = check_packing(A, H)
        o = ilp_opt(A)[0]
        print('N', N, 'n', len(A), '|H|', len(H), 'valid', ok, why, 'pair-maximal', is_pair_maximal(A, H), 'OPT', o, 'excess', o - len(H), flush=True)
        if not ok: continue
        fam = [S for r in range(N + 1) for S in itertools.combinations(range(N), r) if face(A, [xs[i] for i in S], o)]
        maxf = [S for S in fam if not any(set(S) < set(T) for T in fam)]
        print('   maximal faces:', maxf, flush=True)

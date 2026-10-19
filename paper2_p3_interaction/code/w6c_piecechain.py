"""Chain of flip gadgets G* in which consecutive gadgets share a hard piece (K2) of their K4's."""
import sys, itertools, numpy as np
from w6c_lib import *
from w6c_ringK4 import face
GE = [(0,1),(0,4),(0,5),(1,4),(1,5),(4,5),(3,0),(3,1),(6,0),(6,5),(3,6),(3,7)]
GM = [(6,3,7),(4,1,3),(5,0,3),(4,5,6),(6,0,1)]
P1 = (0, 4); P2 = (1, 5)
def build(N, src, dst, flip):
    """gadget g+1's piece `src` (P1 or P2) is identified with gadget g's piece `dst`; flip swaps the two vertices"""
    ident = {}
    sp = P1 if src == 'P1' else P2; dp = P1 if dst == 'P1' else P2
    if flip: dp = dp[::-1]
    for g in range(1, N): 
        for k in range(2): ident[(g, sp[k])] = (g - 1, dp[k])
    idx = {}
    def res(g, v):
        while (g, v) in ident: g, v = ident[(g, v)]
        return (g, v)
    def V(g, v):
        r = res(g, v)
        if r not in idx: idx[r] = len(idx)
        return idx[r]
    members = []; edges = set()
    for g in range(N):
        for (a, b) in GE: edges.add(key(V(g, a), V(g, b)))
        for m in GM: members.append(tuple(V(g, t) for t in m))
    A = adj_from_edges(len(idx), edges)
    return A, members, [V(g, 7) for g in range(N)]
for src in ('P1', 'P2'):
    for dst in ('P1', 'P2'):
        for flip in (False, True):
            for N in (3, 4):
                A, H, xs = build(N, src, dst, flip)
                ok, why = check_packing(A, H)
                if not ok: print(src, dst, flip, N, 'invalid', why); continue
                o = ilp_opt(A)[0]
                if o != len(H): print(src, dst, flip, N, 'notight', o - len(H)); continue
                fam = [S for r in range(N + 1) for S in itertools.combinations(range(N), r) if face(A, [xs[i] for i in S], o)]
                maxf = [S for S in fam if not any(set(S) < set(T) for T in fam)]
                print(src, dst, flip, N, 'n', len(A), 'tight; maximal faces', maxf, flush=True)

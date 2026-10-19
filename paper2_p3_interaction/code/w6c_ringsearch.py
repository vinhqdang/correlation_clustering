"""Search over ways to chain the flip gadget G* into a ring (wire member through a new vertex t, or vertex identification).
Report options for which the ring of N anchors has a lock core of size N (full set not a face, every N-1 subset a face)."""
import sys, itertools, numpy as np, time
from w6c_lib import *
from w6c_ringK4 import face
GV = [0, 1, 3, 4, 5, 6, 7]
GE = [(0,1),(0,4),(0,5),(1,4),(1,5),(4,5),(3,0),(3,1),(6,0),(6,5),(3,6),(3,7)]
GM = [(6,3,7),(4,1,3),(5,0,3),(4,5,6),(6,0,1)]
def build(N, opt):
    idx = {}
    def v(name):
        if name not in idx: idx[name] = len(idx)
        return idx[name]
    ident = {}
    if opt[0] == 'ident':
        u, w = opt[1], opt[2]
    members = []; edges = set()
    def vv(g, t):
        if opt[0] == 'ident' and t == opt[2]:
            return v((opt[1], (g - 1) % N))   # vertex w of gadget g is identified with vertex u of gadget g-1
        return v((t, g))
    for g in range(N):
        for (a, b) in GE: edges.add(key(vv(g, a), vv(g, b)))
        for m in GM: members.append(tuple(vv(g, t) for t in m))
    if opt[0] == 'wire':
        _, u, tp, w = opt
        for g in range(N):
            a = vv(g, u); b = vv((g + 1) % N, w); t = v(('t', g))
            W = [None] * 3; W[tp] = t; rest = [a, b]; k = 0
            for i in range(3):
                if W[i] is None: W[i] = rest[k]; k += 1
            members.append(tuple(W)); edges |= {key(W[0], W[1]), key(W[1], W[2])}
    n = len(idx); A = adj_from_edges(n, edges)
    xs = [idx[(7, g)] for g in range(N)]
    return A, members, xs
def evaluate(N, opt):
    A, H, xs = build(N, opt)
    ok, why = check_packing(A, H)
    if not ok: return ('invalid', why)
    o = ilp_opt(A)[0]
    if o != len(H): return ('notight', o - len(H))
    fam = [S for r in range(N + 1) for S in itertools.combinations(range(N), r) if face(A, [xs[i] for i in S], o)]
    maxf = [S for S in fam if not any(set(S) < set(T) for T in fam)]
    return ('ok', maxf)
if __name__ == '__main__':
    part = int(sys.argv[1]); nparts = int(sys.argv[2])
    opts = [('wire', u, tp, w) for u in (0, 1, 3, 4, 5, 6) for tp in (0, 1, 2) for w in (0, 1, 3, 4, 5, 6)]
    opts += [('ident', u, w) for u in (0, 1, 3, 4, 5, 6) for w in (0, 1, 3, 4, 5, 6)]
    for i, opt in enumerate(opts):
        if i % nparts != part: continue
        r3 = evaluate(3, opt)
        line = f'{opt} N=3 {r3}'
        if r3[0] == 'ok' and r3[1] != [()]:
            r5 = evaluate(5, opt); line += f' | N=5 {r5}'
        print(line, flush=True)

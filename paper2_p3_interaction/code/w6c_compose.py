"""Compose the 5-member flip gadget (w6c_flipgadget.py) into a chain by sharing one vertex between the right interface member of
gadget j and the left interface member of gadget j+1 and a pattern of the 4 cross pairs; test by ILP:
 T1 all locked: pairs (state I1_first, state I2_last) feasible = {(M,M),(R,R)};  T2 one release gives (M,R) ."""
import itertools, sys, numpy as np
from w6c_lib import *
ST = {'L': (0, 1, 0), 'R': (1, 0, 0), 'M': (1, 1, 1)}
G_members = [(6,3,7),(4,1,3),(5,0,3),(4,5,6),(6,0,1)]
G_hard = [(0,4),(1,5)]
X = 7; I1 = 3; I2 = 1       # indices of interface members
def build(N, share, cross_pat):
    """gadgets g=0..N-1 with vertex offset 8*g, vertex v of I2_g (local index share[0] in member) identified with
    vertex of I1_{g+1} (local index share[1]); cross_pat: 4-bit pattern of hard pairs between I2_g and I1_{g+1} other vertices."""
    vid = {}                       # (g, v) -> global
    cnt = 0
    for g in range(N):
        for v in range(8):
            vid[(g, v)] = cnt; cnt += 1
    # identification
    ident = {}
    for g in range(N - 1):
        a = G_members[I2][share[0]]; b = G_members[I1][share[1]]
        ident[(g + 1, b)] = (g, a)
    def res(g, v):
        while (g, v) in ident: g, v = ident[(g, v)]
        return vid[(g, v)]
    newid = {}; nn = 0
    for g in range(N):
        for v in range(8):
            r = res(g, v)
            if r not in newid: newid[r] = nn; nn += 1
    R = lambda g, v: newid[res(g, v)]
    members = []; edges = set()
    for g in range(N):
        for m in G_members:
            t = tuple(R(g, v) for v in m); members.append(t); edges |= {key(t[0], t[1]), key(t[1], t[2])}
        for (a, b) in G_hard: edges.add(key(R(g, a), R(g, b)))
    # cross hard pairs
    for g in range(N - 1):
        o1 = [v for i, v in enumerate(G_members[I2]) if i != share[0]]
        o2 = [v for i, v in enumerate(G_members[I1]) if i != share[1]]
        k = 0
        for a in o1:
            for b in o2:
                if (cross_pat >> k) & 1: edges.add(key(R(g, a), R(g + 1, b)))
                k += 1
    A = adj_from_edges(nn, edges)
    info = dict(members=members, X=[R(g, X) for g in range(N)], I1=[members[g * 5 + I1] for g in range(N)], I2=[members[g * 5 + I2] for g in range(N)])
    return A, info
def feasible(A, H, fixed):
    r = ilp_opt(A, fixed=fixed)
    return r is not None and r[0] == len(H)
def state_fix(mem, s):
    a, b, c = mem; ab, bc, ac = ST[s]
    return {key(a, b): ab, key(b, c): bc, key(a, c): ac}
def locked_fix(n, x):
    return {key(x, t): 0 for t in range(n) if t != x}
def test(N, share, pat, verbose=False):
    A, info = build(N, share, pat); n = len(A); H = info['members']
    ok, why = check_packing(A, H)
    if not ok: return None
    if ilp_opt(A)[0] != len(H): return ('notight',)
    res = {}
    def pairs_for(lockset):
        out = set()
        for s in 'RM':
            for t in 'RM':
                fx = {}
                fx.update(state_fix(info['I1'][0], s)); fx.update(state_fix(info['I2'][-1], t))
                for g in lockset: fx.update(locked_fix(n, info['X'][g]))
                if feasible(A, H, fx): out.add((s, t))
        return out
    allL = pairs_for(range(N))
    rel = [pairs_for([g for g in range(N) if g != j]) for j in range(N)]
    return allL, rel
if __name__ == '__main__':
    N = int(sys.argv[1]); good = []
    for sh0 in range(3):
        for sh1 in range(3):
            for pat in range(16):
                r = test(N, (sh0, sh1), pat)
                if r is None or r[0] == 'notight': continue
                allL, rel = r
                flipok = (('M', 'R') in allL) is False and all(('M', 'R') in x for x in rel) and allL == {('M', 'M'), ('R', 'R')}
                print((sh0, sh1, pat), 'allL', sorted(allL), 'rel', [sorted(x) for x in rel], 'FLIPCHAIN' if flipok else '', flush=True)
                if flipok: good.append((sh0, sh1, pat))
    print('good', good)

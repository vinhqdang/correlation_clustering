"""w11_maj3.py -- K3 scope with the six hard pairs as FLIP VARIABLES (local hard edits allowed): is the relation closed under (M-majority on states, Boolean majority on flips)?"""
import itertools, sys
import numpy as np
from collections import Counter
from w11_maj2 import centre_variants, TOG

def relation_with_flips(members_vs, nv, G_bits):
    verts = range(nv)
    soft = {}
    roles_list = None
    return None

def build_flip_relation(roles, hardpairs, Gbits, nv):
    m = len(roles)
    soft = {}
    for i, (a, b, c) in enumerate(roles):
        for p, k in (((a, b), 0), ((b, c), 1), ((a, c), 2)):
            soft[tuple(sorted(p))] = (i, k)
    hidx = {p: j for j, p in enumerate(hardpairs)}
    states = list(itertools.product((0, 1, 2), repeat=m))
    flips = list(itertools.product((0, 1), repeat=len(hardpairs)))
    R = []
    triples = list(itertools.combinations(range(nv), 3))
    for s in states:
        for f in flips:
            ok = True
            for t in triples:
                cnt = 0
                for p in itertools.combinations(t, 2):
                    if p in soft:
                        i, k = soft[p]; cnt += TOG[s[i]][k]
                    else:
                        j = hidx[p]; cnt += Gbits[j] ^ f[j]
                if cnt == 2: ok = False; break
            if ok: R.append(s + f)
    return R

def Mmaj(x, y, z):
    if x == y or x == z: return x
    if y == z: return y
    return 2

def closed(R, m):
    S = set(R)
    L = list(R)
    for u in L:
        for v in L:
            for w in L:
                o = tuple(Mmaj(u[i], v[i], w[i]) for i in range(m)) + tuple(int(u[i] + v[i] + w[i] >= 2) for i in range(m, len(u)))
                if o not in S: return False
    return True

if __name__ == '__main__':
    A = [0, 1, 3]; B = [0, 2, 4]; C = [1, 2, 5]
    nv = 6
    soft_pairs = set()
    for vs in (A, B, C):
        for p in itertools.combinations(sorted(vs), 2): soft_pairs.add(p)
    hardpairs = [p for p in itertools.combinations(range(nv), 2) if p not in soft_pairs]
    print('hard pairs', hardpairs)
    cnt = Counter(); bad_ex = []
    for roles in itertools.product(*[centre_variants(vs) for vs in (A, B, C)]):
        for Gb in itertools.product((0, 1), repeat=6):
            R = build_flip_relation(roles, hardpairs, Gb, nv)
            if not R: cnt['empty'] += 1; continue
            c = closed(R, 3)
            cnt['closed' if c else 'NOT'] += 1
            if not c and len(bad_ex) < 3: bad_ex.append((roles, Gb, len(R)))
        # only a few role combos to keep time reasonable
        break
    print(cnt); print(bad_ex)

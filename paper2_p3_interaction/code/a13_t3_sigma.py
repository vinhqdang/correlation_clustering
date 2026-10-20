"""a13_t3_sigma.py: structure lemma for Sigma_0(G_N) (all far pairs dropped, ZE states, local hard pairs unedited), by brute force over the
3^(3N) state vectors and by scope analysis of the active triples."""
import itertools, sys
import numpy as np
from a13_lib import *
S2 = {'L': 0, 'R': 1, 'M': 2}
def final_status(inst, A, states, u, v):
    """final edge status of pair uv given ZE states; counted soft pairs edited by state, hard pairs unedited"""
    return None
def sigma0(inst):
    n, mem = inst['n'], inst['members']
    A = adjm(n, inst['E']); cnt = counted_pairs(inst)
    soft_owner = {}
    for i, m in enumerate(mem):
        a, b, c = m
        soft_owner[frozenset((a, b))] = (i, 0); soft_owner[frozenset((b, c))] = (i, 1); soft_owner[frozenset((a, c))] = (i, 2)
    # active triples: all three pairs counted
    act = []
    for tr in itertools.combinations(range(n), 3):
        ps = [frozenset(p) for p in itertools.combinations(tr, 2)]
        if all(p in cnt for p in ps): act.append(tr)
    # constraint scopes
    scopes = {}
    for tr in act:
        ps = [frozenset(p) for p in itertools.combinations(tr, 2)]
        sc = tuple(sorted(set(soft_owner[p][0] for p in ps if p in soft_owner)))
        scopes.setdefault(sc, []).append(tr)
    return act, scopes, soft_owner, A

def enum_sigma(inst, act, soft_owner, A):
    m = len(inst['members'])
    sols = []
    for st in itertools.product(range(3), repeat=m):
        def edge(u, v):
            p = frozenset((u, v))
            if p in soft_owner:
                i, k = soft_owner[p]
                s = A[u][v]
                # state L edits pair 0 (ab), R pair 1 (bc), M pair 2 (ac)
                return 1 - s if st[i] == k else s
            return A[u][v]
        ok = True
        for (x, y, z) in act:
            if edge(x, y) + edge(y, z) + edge(x, z) == 2: ok = False; break
        if ok: sols.append(st)
    return sols

for N in range(1, 5):
    inst = build_GN(N)
    act, scopes, so, A = sigma0(inst)
    nm = len(inst['members'])
    # scope sizes / which members: member index 3(j-1)+{0,1,2} = A_j,B_j,C_j
    blk = lambda i: i//3
    kinds = set()
    for sc in scopes:
        names = tuple(('ABC'[i % 3] + str(blk(i)+1)) for i in sc)
        kinds.add((len(sc), tuple(sorted(set(blk(i) for i in sc))) if len(set(blk(i) for i in sc)) == 1 else 'two', ))
    cross = sorted(set(tuple(('ABC'[i % 3] + str(blk(i)+1)) for i in sc) for sc in scopes if len(set(blk(i) for i in sc)) > 1))
    print('N=%d  #active triples=%d, scopes with members in >1 block: %s' % (N, len(act), cross))
    if nm <= 12:
        sols = enum_sigma(inst, act, so, A)
        # claimed set
        claimed = []
        lab = {'LRL': (0, 1, 0), 'RLR': (1, 0, 1), 'RRM': (1, 1, 2)}
        for ch in itertools.product(lab, repeat=N):
            if all(not (ch[j] == 'RLR' and ch[j+1] == 'LRL') for j in range(N-1)):
                claimed.append(tuple(x for c in ch for x in lab[c]))
        print('   |Sigma_0| brute=%d claimed=%d equal=%s' % (len(sols), len(claimed), set(sols) == set(claimed)), flush=True)

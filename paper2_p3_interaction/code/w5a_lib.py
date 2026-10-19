"""Wave 5 prover A: checks of the pointer-anatomy lemma (T1/T2) and the tiny-tiny locality corollary.
Reuses a4_lib (graph = 0/1 matrix, packed = soft).  No git.
anatomy(A, pk, lab) -> dict of violations / counts for one partition lab (cluster labels)."""
import itertools
import numpy as np
from a4_lib import *

def anatomy(A, pk, lab):
    n = len(A)
    pid, U, F, I = pieces_of(A, pk, lambda i, j: lab[i] == lab[j])
    u = len(U); Iset = set(I)
    res = dict(viol=[], t1=0, t2=0, tiny_pairs=0, u=u)
    cl = {}
    for v in range(n): cl.setdefault(lab[v], []).append(v)
    # L1/T1/T2 for every inserted hard pair
    for (a, b) in I:
        C = cl[lab[a]]
        zs = [z for z in C if z not in (a, b) and A[a][z] and A[z][b]]
        if not zs: res['viol'].append(('diam', a, b)); continue
        soft2 = [z for z in zs if pk[a, z] and pk[z, b]]
        if soft2:
            res['t1'] += 1
            # the two packed pairs lie in different members (else ab packed) -- automatic since ab is hard
        else:
            ok = False
            for z in zs:
                if (not pk[a, z] and pid[a] == pid[z] and pk[z, b]) or (not pk[z, b] and pid[b] == pid[z] and pk[a, z]): ok = True
            if ok: res['t2'] += 1
            else: res['viol'].append(('allhard', a, b))
    # tiny-tiny corollary: pieces in one cluster with no soft pair between them
    for c, C in cl.items():
        ps = {}
        for v in C: ps.setdefault(pid[v], []).append(v)
        for P, Q in itertools.combinations(list(ps.values()), 2):
            if any(pk[a, b] for a in P for b in Q): continue
            res['tiny_pairs'] += 1
            PQ = set(P) | set(Q)
            for a in P:
                for b in Q:
                    if key(a, b) not in Iset: res['viol'].append(('cross-not-inserted', a, b))
                    if not any(A[a][z] and A[z][b] and pk[a, z] and pk[z, b] for z in C if z not in PQ):
                        res['viol'].append(('tiny-no-soft-z', a, b))
    # near-twin / agreement lemma: z,b same cluster => |N_h(z) \ (N[b] u Q(b))| <= |U|
    for c, C in cl.items():
        for z, b in itertools.permutations(C, 2):
            M = [a for a in range(n) if a not in (z, b) and not pk[a, z] and not pk[a, b] and A[a][z] != A[a][b]]
            if len(M) > u: res['viol'].append(('M-bound', z, b, len(M), u))
    return res

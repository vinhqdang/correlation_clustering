"""Random test of Theorem A + Cor A1 with independent brute force (all set partitions), n<=10.
usage: python3 a10_closure.py seed count"""
import sys, random, itertools
import numpy as np
from a10_lib import *

def apply(A, H, assign):
    B = A.copy()
    for i, s in assign.items():
        a, b, c = H[i]; x, y = [(a, b), (b, c), (a, c)][s]
        B[x, y] = B[y, x] = not B[x, y]
    return B

def pc_feasible(A, H):
    n = len(A); m = len(H); o = owner_map(H)
    D = [set(range(3)) for _ in range(m)]
    R = {(i, j): np.ones((3, 3), bool) for i in range(m) for j in range(m) if i != j}
    for t in itertools.combinations(range(n), 3):
        ow = sorted({o[p] for p in map(frozenset, itertools.combinations(t, 2)) if p in o})
        if len(ow) > 2: continue     # soft triangle: ternary, ignored by PC
        ok = {}
        for st in itertools.product(range(3), repeat=len(ow)):
            B = apply(A, H, dict(zip(ow, st)))
            x, y, z = t
            ok[st] = int(B[x, y]) + int(B[y, z]) + int(B[x, z]) != 2
        if len(ow) == 0:
            if not ok[()]: return False
        elif len(ow) == 1:
            D[ow[0]] &= {s for (s,), v in ok.items() if v}
        else:
            i, j = ow
            for (s, u), v in ok.items():
                if not v: R[(i, j)][s, u] = False; R[(j, i)][u, s] = False
    ch = True
    while ch:
        ch = False
        for (i, j), M in R.items():
            for s in range(3):
                if s not in D[i]: M[s, :] = False
                if s not in D[j]: M[:, s] = False
            for k in range(m):
                if k in (i, j): continue
                comp = (R[(i, k)].astype(int) @ R[(k, j)].astype(int)) > 0
                new = M & comp
                if (new != M).any(): M[:] = new; ch = True
            for s in range(3):
                if s in D[i] and not M[s].any(): D[i].discard(s); ch = True
        if any(len(d) == 0 for d in D): return False
        for (i, j), M in R.items():
            pass
    return all(len(d) > 0 for d in D) and all(M.any() for M in R.values())

def main():
    seed, count = int(sys.argv[1]), int(sys.argv[2])
    rng = random.Random(seed)
    st = dict(tot=0, stf=0, stf_closed=0, stf_2dec=0, stf_nonempty=0, stf_big=0, pc_agree=0, pc_dis=0,
              sf=0, sf_notclosed=0, sf_not2dec=0, cmaxhist={}, pc_sf_wrong=0, sf_empty_ok=0)
    for it in range(count):
        n = rng.randint(5, 10); m = rng.randint(2, 6)
        g = planted(n, m, rng, big=rng.random() < 0.6)
        if g is None: continue
        A, H = g
        # noise: flip some random hard pairs
        if rng.random() < 0.5:
            o = owner_map(H)
            for _ in range(rng.randint(1, 2)):
                x, y = rng.sample(range(n), 2)
                if frozenset((x, y)) not in o: A[x, y] = A[y, x] = not A[x, y]
        if not valid(A, H): continue
        if len(partitions(n)) > 120000: continue
        sig, mn = sigma_bruteforce(A, H)
        # NB: with noise the min partition cost can exceed |H|; sig then empty
        st['tot'] += 1
        c = cmax(n, H); st['cmaxhist'][c] = st['cmaxhist'].get(c, 0) + 1
        stri = soft_triangles(n, H)
        feas = len(sig) > 0
        if not stri:
            st['stf'] += 1
            if feas: st['stf_nonempty'] += 1
            if len(sig) >= 3: st['stf_big'] += 1
            if closed_all_rules(sig): st['stf_closed'] += 1
            else: print("COUNTEREXAMPLE closure", n, A.astype(int).tolist(), H, sorted(sig))
            if two_decomposable(sig, len(H)): st['stf_2dec'] += 1
            else: print("COUNTEREXAMPLE 2dec", H)
            if pc_feasible(A, H) == feas: st['pc_agree'] += 1
            else: st['pc_dis'] += 1; print("PC DISAGREE", n, H, feas)
        else:
            st['sf'] += 1
            if not closed_all_rules(sig): st['sf_notclosed'] += 1
            if not two_decomposable(sig, len(H)):
                st['sf_not2dec'] += 1; st.setdefault('not2dec_by_c',{}); st['not2dec_by_c'][c]=st['not2dec_by_c'].get(c,0)+1
            if not closed_all_rules(sig): st.setdefault('notclosed_by_c',{}); st['notclosed_by_c'][c]=st['notclosed_by_c'].get(c,0)+1
            if pc_feasible(A, H) != feas: st['pc_sf_wrong'] += 1
    print(st)
if __name__=="__main__": main()

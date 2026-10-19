"""Audit 7: brute force of W1, W2 (+tightness), Cor LB, F1, DEL, BAD on all partitions of cost <= |H|+3."""
import sys, random, math, itertools
import numpy as np
import a7_lib as L

def run(seed, ninst, nrange, slack=3):
    rng = random.Random(seed)
    st = dict(inst=0, parts=0, w1=0, w1_viol=0, w2_viol=0, lb_viol=0, f1_viol=0, f1_checks=0, del_viol=0, del_checks=0,
              bad_viol=0, bad_checks=0, maxratio=0.0, ratio1=0, nonmax=0)
    ex = None
    for it in range(ninst):
        n = rng.choice(nrange)
        A, H = L.gen_instance(n, rng)
        if rng.random() < .25 and len(H) > 1:   # non-maximal packing: W1/CP do not need maximality
            H = H[:-1]; st['nonmax'] += 1
        if not H: continue
        S = L.soft_of(n, H); Dl, Gm = L.delta_gamma(A, S)
        c = max(sum(1 for t in H if v in t) for v in range(n)); D = 2 * c
        P = L.partitions(n); cs, X, pl = L.costs(A, P)
        idx = {p: k for k, p in enumerate(pl)}
        E = np.array([bool(A[i][j]) for i, j in pl])
        ok = np.where(cs <= len(H) + slack)[0]
        st['inst'] += 1
        # static: hard components for DEL; N_h
        Nh = [set(y for y in range(n) if A[v][y] and not S[v, y]) for v in range(n)]
        for r in ok:
            st['parts'] += 1
            row = X[r]; ed = row != E
            cost = int(cs[r]); eps = cost - len(H)
            U = set(pl[k] for k in np.where(ed)[0] if not S[pl[k]])
            u = len(U)
            assert cost >= len(H)
            # W1 + W2 weights
            Wt = 0; wt = {}
            for z, b in L.prs(n):
                if not S[z, b]: continue
                merged = bool(row[idx[z, b]])
                ys = Dl[z, b] if merged else Gm[z, b]
                for y in ys:
                    p1 = tuple(sorted((z, y))); p2 = tuple(sorted((b, y)))
                    st['w1'] += 1
                    if p1 not in U and p2 not in U: st['w1_viol'] += 1; ex = ('W1', A.copy(), H, r)
                wt[z, b] = len(ys); Wt += len(ys)
            if Wt > 2 * D * u: st['w2_viol'] += 1; ex = ('W2', A.copy(), H, r)
            if u > 0: 
                rt = Wt / (2 * D * u); st['maxratio'] = max(st['maxratio'], rt)
                if abs(rt - 1) < 1e-12: st['ratio1'] += 1
            # c_i, CorLB
            ci = []
            for t in H:
                a, b_, c_ = t
                ci.append(sum(1 for p in p3p_list(t) if ed[idx[p]]))
            if min(ci) < 1: raise SystemExit("c_i<1")
            assert cost == sum(ci) + u, "CP(b)"
            if cost < sum(ci) + math.ceil(Wt / (2 * D)): st['lb_viol'] += 1
            # F1 (needs |S| <= |H|+l : take l = eps, and l = eps+1,2)
            for l in range(eps, eps + 3):
                for (z, b), w in wt.items():
                    merged = bool(row[idx[z, b]])
                    if merged:
                        st['f1_checks'] += 1
                        if len(Dl[z, b]) > l: st['f1_viol'] += 1
                    else:
                        st['f1_checks'] += 1
                        if len(Gm[z, b]) > l: st['f1_viol'] += 1
            # DEL
            for (p, q) in U:
                if A[p][q]:   # deleted hard edge
                    st['del_checks'] += 1
                    if u < 1 + len(Nh[p] & Nh[q]): st['del_viol'] += 1
            # BAD
            nbad = 0
            for t in H:
                x, y, z = t
                states = [ ({(y, z)}), ({(x, y)}), ({(x, y), (y, z), (x, z)}) ]   # merged soft pairs in states: del xy, del yz, ins xz
                states = [{tuple(sorted(p)) for p in s} for s in states]
                allsoft = [tuple(sorted(p)) for p in ((x, y), (y, z), (x, z))]
                bad = True
                for s in states:
                    w = 0
                    for pr in allsoft:
                        w += len(Dl[pr]) if pr in s else len(Gm[pr])
                    if w == 0: bad = False
                if bad: nbad += 1
            st['bad_checks'] += 1
            if nbad > 2 * c * eps: st['bad_viol'] += 1; ex = ('BAD', A.copy(), H, r)
    return st, ex

def p3p_list(t):
    a, b, c = t
    return [tuple(sorted(p)) for p in ((a, b), (b, c), (a, c))]

if __name__ == "__main__":
    seed = int(sys.argv[1]); ninst = int(sys.argv[2])
    st, ex = run(seed, ninst, [5, 6, 6, 7, 7, 8, 8, 9] if len(sys.argv) < 4 else [int(x) for x in sys.argv[3].split(',')])
    print(st)
    if ex: print("EXAMPLE VIOLATION", ex[0], ex[2], ex[1].tolist())

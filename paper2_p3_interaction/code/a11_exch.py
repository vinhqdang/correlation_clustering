"""Audit 11: exchange calculus, twin lemma, witness dichotomy, displacement, hard components of uncovered vertices.
Ground truth: all optima by full partition enumeration.  usage: python3 a11_exch.py SEED TRIALS"""
import sys, random, itertools
from a11_lib import *

def clusters_of(g, S):
    A = [[g.e(a, b) != (K(a, b) in S) if a != b else False for b in range(g.n)] for a in range(g.n)]
    seen = {}; cl = []
    for v in range(g.n):
        if v in seen: continue
        comp = [v]; seen[v] = 1; st = [v]
        while st:
            x = st.pop()
            for y in range(g.n):
                if A[x][y] and y not in seen: seen[y] = 1; comp.append(y); st.append(y)
        cl.append(set(comp))
    return cl

def main(seed, trials):
    rng = random.Random(seed)
    st = dict(inst=0, optima=0, hardedits=0, twin_viol=0, twin_pairs=0, noP3=0, witness_soft0=0, dich_a=0, dich_b=0, dich_viol=0,
              exch_bad=0, exch_tw_bad=0, exch_notpm=0, exch_pm=0, mid_notmiddle=0, both_unc=0, l25=0, l25_bad=0, l25_out=0, l25_in=0,
              l26=0, l26_bad=0, dich_a_notlocal=0)
    for _ in range(trials):
        n = rng.randint(6, 10); p = rng.choice([0.25, 0.35, 0.5])
        g = Gr(n, [(a, b) for a in range(n) for b in range(a + 1, n) if rng.random() < p])
        if not g.p3(): continue
        H = random_packing(g, rng)
        assert is_packing(g, H) and pair_maximal(g, H)
        o, opts = opt_and_optima(g)
        st['inst'] += 1
        own = soft_owner(H)
        cov = {v for t in H for v in t}
        I = I_graph(H); tw0 = tw_exact(I)
        P3 = g.p3()
        # twins
        for a in range(n):
            for b in range(a + 1, n):
                if g.adj[a] | {a} == g.adj[b] | {b}:
                    st['twin_pairs'] += 1
                    if any((a, b) in S for S in opts): st['twin_viol'] += 1
        # Lemma 2.6 (outside vertices): hard component of r
        hardE = [pq for pq in g.E if pq not in own]
        for r in range(n):
            if r in cov: continue
            st['l26'] += 1
            comp = {r}; stack = [r]
            while stack:
                x = stack.pop()
                for (u, v) in hardE:
                    if u == x and v not in comp: comp.add(v); stack.append(v)
                    if v == x and u not in comp: comp.add(u); stack.append(u)
            ok = comp == {r} | g.adj[r]
            for u, v in itertools.combinations(sorted(comp), 2):
                if not (K(u, v) in g.E and K(u, v) not in own) and K(u, v) not in own: ok = False
            if not ok: st['l26_bad'] += 1
        for S in opts:
            st['optima'] += 1
            U = hard_of(H, S)
            cl = clusters_of(g, S)
            clof = {v: i for i, C in enumerate(cl) for v in C}
            for e in U:
                st['hardedits'] += 1
                a, b = e
                if a not in cov and b not in cov: st['both_unc'] += 1
                wit = [t for t in P3 if e in pr(t)]
                if not wit: st['noP3'] += 1; continue
                nsoft = [sum(1 for q in pr(t) if q in own) for t in wit]
                if any(k == 0 or k == 3 for k in nsoft): st['dich_viol'] += 1
                two = [t for t, k in zip(wit, nsoft) if k == 2]
                if two:
                    st['dich_a'] += 1
                    for t in two:
                        sp = [q for q in pr(t) if q in own]
                        P1, P2 = own[sp[0]], own[sp[1]]
                        z = set(sp[0]) & set(sp[1])
                        z = next(iter(z))
                        ok = P1 != P2 and z in H[P1] and z in H[P2] and a in (set(H[P1]) | set(H[P2])) and b in (set(H[P1]) | set(H[P2]))
                        ok = ok and ((a in H[P1] and b in H[P2]) or (a in H[P2] and b in H[P1]))
                        if not ok: st['dich_viol'] += 1
                        if z != t[1]:
                            st['mid_notmiddle'] += 1
                            if e not in g.E: st['mid_notmiddle_ins'] = st.get('mid_notmiddle_ins', 0) + 1
                else:
                    st['dich_b'] += 1
                    # every witness exactly one soft pair; exchange
                    for t in wit:
                        s = [q for q in pr(t) if q in own][0]
                        Q = own[s]
                        H2 = [x for i, x in enumerate(H) if i != Q] + [t]
                        okp = is_packing(g, H2) and len(H2) == len(H)
                        if not okp: st['exch_bad'] += 1; continue
                        o2 = soft_owner(H2)
                        if e not in o2: st['exch_bad'] += 1
                        tw2 = tw_exact(I_graph(H2))
                        if tw2 > tw0 + 1: st['exch_tw_bad'] += 1
                        if pair_maximal(g, H2): st['exch_pm'] += 1
                        else: st['exch_notpm'] += 1
                        # Lemma 2.5 (inserted e, z kept in cluster)
                        if e not in g.E:
                            z = t[1]
                            if clof[z] == clof[a] == clof[b] and z != a and z != b:
                                st['l25'] += 1
                                Qs = pr(H[Q]); edits_in_Q = [q for q in Qs if q in S]
                                # soft step is s (contains z and a or b)
                                if len(edits_in_Q) != 1: st['l25_bad'] += 1; continue
                                f = edits_in_Q[0]
                                newU = set(hard_of(H2, S))
                                if newU != (set(U) - {e}) | {f}: st['l25_bad'] += 1; continue
                                r_ = [v for v in H[Q] if v not in s][0]
                                if f[0] != r_ and f[1] != r_: st['l25_bad'] += 1
                                cov2 = {v for tt in H2 for v in tt}
                                if r_ in cov2: st['l25_in'] += 1
                                else: st['l25_out'] += 1
    print(seed, st)

if __name__ == '__main__':
    main(int(sys.argv[1]), int(sys.argv[2]))

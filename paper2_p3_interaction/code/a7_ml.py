"""Audit 7: Lemma D, Rule 6.2 on tight components, Theorem 3 counting (|A|<=2l', |B|<=l-l', r<=l+l').  Own code."""
import sys, random, itertools
import numpy as np, networkx as nx
import a7_lib as L

def comps_of(H):
    I = nx.Graph(); I.add_nodes_from(range(len(H)))
    for i, j in itertools.combinations(range(len(H)), 2):
        if set(H[i]) & set(H[j]): I.add_edge(i, j)
    return [sorted(c) for c in nx.connected_components(I)]

def Wset(H, comp): return sorted(set(v for i in comp for v in H[i]))

def is_cluster(A, Vs):
    for a, b, c in itertools.combinations(Vs, 3):
        if A[a][b] + A[b][c] + A[a][c] == 2: return False
    return True

def clusters_of(Aw, Vs):
    G = nx.Graph(); G.add_nodes_from(Vs); G.add_edges_from((i, j) for i, j in itertools.combinations(Vs, 2) if Aw[i][j])
    return [set(c) for c in nx.connected_components(G)]

def R_ok(A, Aw, Vs, W):
    """R-condition for solution graph Aw on W: every outside y has N(y)&W empty or exactly one cluster of Aw[W]."""
    cl = clusters_of(Aw, W)
    Wset_ = set(W)
    for y in range(len(A)):
        if y in Wset_: continue
        N = {u for u in W if A[y][u]}
        if N and N not in cl: return False
    return True

def tight_solutions(A, H, comp):
    """all S (set of pairs) with exactly one edited pair per member of comp and G[W] xor S a cluster graph."""
    W = Wset(H, comp); res = []
    pl = []
    for i in comp:
        x, y, z = H[i]; pl.append([tuple(sorted((x, y))), tuple(sorted((y, z))), tuple(sorted((x, z)))])
    for ch in itertools.product(range(3), repeat=len(comp)):
        S = [pl[k][ch[k]] for k in range(len(comp))]
        if len(set(S)) < len(S): continue
        B = A.copy()
        for a, b in S: B[a][b] ^= 1; B[b][a] ^= 1
        if is_cluster(B, W): res.append((S, B))
    return W, res

def apply_S(A, S):
    B = A.copy()
    for a, b in S: B[a][b] ^= 1; B[b][a] ^= 1
    return B

def gen(n, rng):
    while True:
        A = L.gen_graph(n, rng)
        ts = L.ind_p3s(A)
        if not ts: continue
        mode = rng.choice(["disjoint", "pair", "pair"])
        rng.shuffle(ts); H = []; used = set(); usedv = set()
        for t in ts:
            ps = L.p3p(t)
            if ps & used: continue
            if mode == "disjoint" and set(t) & usedv: continue
            H.append(t); used |= ps; usedv |= set(t)
            if rng.random() < .15: break      # sometimes non-maximal packing (Thm 3 holds for any packing)
        if H: return A, H

def main(seed, N, nr):
    rng = random.Random(seed)
    st = dict(inst=0, apps=0, sound_fail=0, pipeline=0, opt_drop_fail=0, A_fail=0, B_fail=0, r_fail=0, notfullyreduced=0, tightcnt=0, irredcnt=0,
              nontight=0, tight_vs_tau_fail=0, lemD_fail=0, lemD_checks=0, lemD_viol_needed=0, multi=0, maxA=0)
    ex = []
    for it in range(N):
        n = rng.choice(nr); A, H = gen(n, rng)
        opt = L.opt_subset(A); cm = comps_of(H)
        info = []
        for comp in cm:
            W, sols = tight_solutions(A, H, comp)
            tau = L.opt_subset(A[np.ix_(W, W)])
            p = len(comp); tight = (tau == p)
            if tight != (len(sols) > 0): st['tight_vs_tau_fail'] += 1; ex.append(('tight', A.tolist(), H))
            good = [(S, B) for S, B in sols if R_ok(A, B, None, W)] if tight else []
            info.append(dict(comp=comp, W=W, tau=tau, p=p, tight=tight, sols=sols, good=good))
            st['tightcnt'] += tight
        st['inst'] += 1
        if len(cm) > 1: st['multi'] += 1
        # soundness of Rule 6.2 for EVERY R-satisfying tight S_F of every component
        for d in info:
            for S, B in d['good']:
                st['apps'] += 1
                if L.opt_subset(B) + len(S) != opt: st['sound_fail'] += 1; ex.append(('sound', A.tolist(), H, S))
        # pipeline: reduce all reducible tight comps (one S each), in order
        Acur = A.copy(); remaining = []
        for d in info:
            if d['tight'] and d['good']:
                S, B = d['good'][0]; Acur = apply_S(Acur, S)
            else: remaining.append(d)
        sump = sum(d['p'] for d in info if d['tight'] and d['good'])
        optr = L.opt_subset(Acur)
        if optr != opt - sump: st['opt_drop_fail'] += 1; ex.append(('optdrop', A.tolist(), H))
        # remaining comps still irreducible in the reduced graph (Lemma K2)
        Aset = []
        for d in remaining:
            if d['tight']:
                W, sols = tight_solutions(Acur, H, d['comp'])
                good = [1 for S, B in sols if R_ok(Acur, B, None, W)]
                if good: st['notfullyreduced'] += 1; ex.append(('notred', A.tolist(), H))
        Acnt = sum(1 for d in remaining if d['tight']); Bcnt = sum(1 for d in remaining if not d['tight'])
        st['irredcnt'] += Acnt; st['nontight'] += Bcnt; st['maxA'] = max(st['maxA'], Acnt)
        for l in range(max(opt - len(H), 0), max(opt - len(H), 0) + 3):
            lp = l - sum(d['tau'] - d['p'] for d in info)
            st['pipeline'] += 1
            if lp < 0: st['A_fail'] += 1; continue      # cannot happen for YES instance
            if Acnt > 2 * lp: st['A_fail'] += 1; ex.append(('A', A.tolist(), H, l))
            if Bcnt > l - lp: st['B_fail'] += 1
            if Acnt + Bcnt > l + lp: st['r_fail'] += 1
        # Lemma D : for ALL partitions (n<=8) or random ones
        if n <= 8 and rng.random() < .5:
            P = L.partitions(n)
            for d in info[:2]:
                W = d['W']; Wsets = set(W)
                for row in P[:: max(1, len(P) // 400)]:
                    Sset = set()
                    for i, j in L.prs(n):
                        if (row[i] == row[j]) != bool(A[i][j]): Sset.add((i, j))
                    B = apply_S(A, Sset)
                    for y in range(n):
                        if y in Wsets: continue
                        N = {u for u in W if A[y][u]}
                        # S_F = S restricted to W
                        BW = A.copy()
                        for a, b in Sset:
                            if a in Wsets and b in Wsets: BW[a][b] ^= 1; BW[b][a] ^= 1
                        cl = clusters_of(BW, W)
                        viol = bool(N) and (N not in cl)
                        has = any((min(y, u), max(y, u)) in Sset for u in W)
                        st['lemD_checks'] += 1; st['lemD_viol_needed'] += viol
                        if viol and not has: st['lemD_fail'] += 1; ex.append(('lemD', A.tolist(), H))
    return st, ex

if __name__ == "__main__":
    s, ex = main(int(sys.argv[1]), int(sys.argv[2]), [int(x) for x in sys.argv[3].split(',')])
    print(s)
    for e in ex[:3]: print(e)

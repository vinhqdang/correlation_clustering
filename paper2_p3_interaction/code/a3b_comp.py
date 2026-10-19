# Wave-3 audit of Theorem K: component-level packing lower bound for Cluster Editing.
# Exhaustive partitions (vectorised) for n <= 10.  Run: python3 -I a3b_comp.py [n_inst]
import itertools, random, sys
import numpy as np

def rgs(n):
    out = []
    def rec(i, cur, m):
        if i == n:
            out.append(tuple(cur)); return
        for b in range(m + 1):
            cur.append(b); rec(i + 1, cur, max(m, b + 1) if b == m else m); cur.pop()
    rec(0, [], 0)
    return np.array(out, dtype=np.int8)

PARTS = {}
def parts(n):
    if n not in PARTS: PARTS[n] = rgs(n)
    return PARTS[n]

def pair_list(n): return list(itertools.combinations(range(n), 2))

def edit_matrix(n, adj, P):
    """N x npairs bool: partition P edits pair (i,j)."""
    pl = pair_list(n)
    same = np.stack([P[:, i] == P[:, j] for i, j in pl], 1)
    e = np.array([adj[i][j] for i, j in pl], dtype=bool)
    return same != e[None, :], pl

def p3s(n, adj):
    out = []
    for b in range(n):
        for a in range(n):
            for c in range(a + 1, n):
                if a != b and c != b and adj[a][b] and adj[b][c] and not adj[a][c]:
                    out.append((a, b, c))
    return out

def packing(n, adj, rng, keep=1.0):
    L = p3s(n, adj); rng.shuffle(L); used = set(); H = []
    for (a, b, c) in L:
        ps = {frozenset(x) for x in ((a, b), (b, c), (a, c))}
        if ps & used or rng.random() > keep: continue
        used |= ps; H.append((a, b, c))
    return H

def components(H):
    comp = []  # list of lists of P3 indices; adjacency = share a vertex
    seen = [False] * len(H)
    for i in range(len(H)):
        if seen[i]: continue
        st = [i]; seen[i] = True; cur = []
        while st:
            x = st.pop(); cur.append(x)
            for y in range(len(H)):
                if not seen[y] and set(H[x]) & set(H[y]): seen[y] = True; st.append(y)
        comp.append(cur)
    return comp

def opt_all(n, adj):
    P = parts(n)
    E, pl = edit_matrix(n, adj, P)
    cost = E.sum(1)
    return P, E, pl, cost

def sub_adj(W, adj): return [[adj[a][b] for b in W] for a in W]

def tau_and_opts(W, adj):
    W = sorted(W); m = len(W)
    P, E, pl, cost = opt_all(m, sub_adj(W, adj))
    t = int(cost.min())
    return t, [(P[i], pl) for i in np.nonzero(cost == t)[0]], W

def reducible_S(W, adj, n):
    """Return (tau, S as set of global pairs) for first optimal S of G[W] meeting Rule 6.2 conditions, else (tau, None)."""
    t, opts, Ws = tau_and_opts(W, adj)
    X = [v for v in range(n) if v not in set(Ws)]
    for lab, pl in opts:
        ok = True
        for ii in range(len(Ws)):
            for jj in range(ii + 1, len(Ws)):
                u, v = Ws[ii], Ws[jj]
                Nu = {x for x in X if adj[u][x]}; Nv = {x for x in X if adj[v][x]}
                if lab[ii] == lab[jj]:
                    if Nu != Nv: ok = False; break
                else:
                    if Nu & Nv: ok = False; break
            if not ok: break
        if ok:
            S = set()
            for (i, j) in pl:
                if (lab[i] == lab[j]) != bool(adj[Ws[i]][Ws[j]]): S.add((Ws[i], Ws[j]))
            return t, S
    return t, None

def run(ninst, seed):
    rng = random.Random(seed)
    st = dict(inst=0, A_opt_ge=0, A_lprime_le=0, A_tau_ge=0, B_ident=0, B_inside_ge=0,
              C_applied=0, C_drop_exact=0, C_bound=0, C_nontrivial=0, C_ratio2=0, D_lemma32=0,
              neg_lprime=0, multi_vertex_share=0, viol=[])
    for it in range(ninst):
        n = rng.choice([7, 8, 9, 10]); p = rng.choice([0.2, 0.3, 0.4, 0.5, 0.6])
        adj = [[0] * n for _ in range(n)]
        for i, j in pair_list(n):
            if rng.random() < p: adj[i][j] = adj[j][i] = 1
        H = packing(n, adj, rng, keep=rng.choice([1.0, 1.0, 0.6]))
        if not H: continue
        st['inst'] += 1
        comps = components(H)
        Ws = [sorted({v for i in c for v in H[i]}) for c in comps]
        if any(len(c) > 1 for c in comps): st['multi_vertex_share'] += 1
        P, E, pl, cost = opt_all(n, adj)
        OPT = int(cost.min()); idx = {q: k for k, q in enumerate(pl)}
        taus = [tau_and_opts(W, adj)[0] for W in Ws]
        # (1) lower bound facts
        for c, t in zip(comps, taus):
            if t >= len(c): st['A_tau_ge'] += 1
            else: st['viol'].append(('tau<|C|', it))
        if OPT >= sum(taus): st['A_opt_ge'] += 1
        else: st['viol'].append(('OPT<sum tau', it))
        l = OPT - len(H); lp = OPT - sum(taus)
        if lp <= l: st['A_lprime_le'] += 1
        else: st['viol'].append(('lp>l', it))
        # (2) decomposition for EVERY partition, hence every optimal one
        inside = np.stack([E[:, [idx[q] for q in itertools.combinations(W, 2)]].sum(1) for W in Ws], 1)
        inW = np.zeros(len(pl), dtype=bool)
        for W in Ws:
            for q in itertools.combinations(W, 2): inW[idx[q]] = True
        ext = E[:, ~inW].sum(1)
        if (inside.sum(1) + ext == cost).all(): st['B_ident'] += 1
        else: st['viol'].append(('ident', it))
        if (inside >= np.array(taus)[None, :]).all(): st['B_inside_ge'] += 1
        else: st['viol'].append(('inside<tau', it))
        # (3) reduction at component level
        G = [r[:] for r in adj]; cur = list(range(len(Ws))); curOPT = OPT; changed = True
        while changed:
            changed = False
            for j in list(cur):
                t, S = reducible_S(Ws[j], G, n)
                if S is not None:
                    for (a, b) in S: G[a][b] ^= 1; G[b][a] ^= 1
                    _, _, _, c2 = opt_all(n, G); newOPT = int(c2.min())
                    st['C_applied'] += 1
                    if newOPT == curOPT - t: st['C_drop_exact'] += 1
                    else: st['viol'].append(('drop', it))
                    curOPT = newOPT; cur.remove(j); changed = True
        if cur:
            tt = [tau_and_opts(Ws[j], G)[0] for j in cur]
            lp2 = curOPT - sum(tt)
            if lp2 < 0: st['neg_lprime'] += 1
            if len(cur) <= 2 * lp2: st['C_bound'] += 1
            else: st['viol'].append(('count>2lp', it, len(cur), lp2))
            st['C_nontrivial'] += 1
            if len(cur) == 2 * lp2: st['C_ratio2'] += 1
            # Lemma 3.2 property on every optimal partition of reduced G
            P2, E2, pl2, c2 = opt_all(n, G); idx2 = {q: k for k, q in enumerate(pl2)}
            opt_rows = np.nonzero(c2 == curOPT)[0]
            ok = True
            for r in opt_rows:
                pa = pb = 0
                for j, t in zip(cur, tt):
                    ins = sum(E2[r, idx2[q]] for q in itertools.combinations(Ws[j], 2))
                    exts = any(E2[r, idx2[q]] for q in pl2 if (q[0] in Ws[j]) != (q[1] in Ws[j]))
                    if ins > t: pa += 1
                    elif exts: pb += 1
                    else: ok = False
                if curOPT < sum(tt) + pa + pb / 2: ok = False
            if ok: st['D_lemma32'] += 1
            else: st['viol'].append(('lemma3.2', it))
    return st

if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    st = run(N, 2024)
    for k, v in st.items(): print(k, v if k != 'viol' else (len(v), v[:5]))

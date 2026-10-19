"""Audit 7: Theorem P (pruned omega theorem).  Brute force + own tree-decomposition DP.
For random pair-maximal (G,H), l=0..3:
  Gamma^l = hard edges in V' + soft pairs in V' with |Delta|<=l ;  V' = V - BIG ;
  f~(C) = C(|C|,2) - 2 e^+(C) , e^+ = G-edges of C that are edges of Gamma^l.
  m_all, m_conn by brute force over partitions of V' (cluster size <= Z); m_dp by DP.
  Criterion  |E| + sum_{K in BIG comps} f(K) + m  <= |H|+l  <=>  OPT <= |H|+l, tested for m_all, m_conn, m_dp."""
import sys, random, itertools
import numpy as np, networkx as nx
from networkx.algorithms.approximation import treewidth_min_fill_in
import a7_lib as L

def setup(A, H, l):
    n = len(A); S = L.soft_of(n, H)
    c = max(sum(1 for t in H if v in t) for v in range(n)); D = 2 * c
    B = 2 * D + 3 * l; Z = B + D + 1
    hard = (A.astype(bool)) & ~S
    dh = hard.sum(1)
    # hard components
    G = nx.Graph(); G.add_nodes_from(range(n)); G.add_edges_from((i, j) for i, j in L.prs(n) if hard[i, j])
    comps = list(nx.connected_components(G))
    BIG = set()
    for K in comps:
        if any(dh[v] > B for v in K): BIG |= K
    Vp = [v for v in range(n) if v not in BIG]
    Dl, Gm = L.delta_gamma(A, S)
    GL = nx.Graph(); GL.add_nodes_from(Vp)
    for i, j in L.prs(n):
        if i in BIG or j in BIG: continue
        if hard[i, j] or (S[i, j] and len(Dl[i, j]) <= l): GL.add_edge(i, j)
    return dict(n=n, S=S, D=D, B=B, Z=Z, BIG=BIG, Vp=Vp, GL=GL, comps=comps, Dl=Dl)

def f_true(A, C):
    C = list(C); s = len(C); e = sum(A[i][j] for i, j in itertools.combinations(C, 2))
    return s * (s - 1) // 2 - 2 * e

def f_tilde(A, GL, C):
    C = list(C); s = len(C); e = sum(1 for i, j in itertools.combinations(C, 2) if A[i][j] and GL.has_edge(i, j))
    return s * (s - 1) // 2 - 2 * e

def brute_m(A, st):
    Vp = st['Vp']; GL = st['GL']; Z = st['Z']; m = len(Vp)
    if m == 0: return 0, 0
    P = L.partitions(m); best_all = 10 ** 9; best_conn = 10 ** 9
    cache = {}
    for row in P:
        blocks = {}
        for idx, lab in enumerate(row): blocks.setdefault(int(lab), []).append(Vp[idx])
        if any(len(b) > Z for b in blocks.values()): continue
        tot = 0; conn = True
        for b in blocks.values():
            key = frozenset(b)
            if key not in cache: cache[key] = (f_tilde(A, GL, b), nx.is_connected(GL.subgraph(b)))
            fv, cv = cache[key]; tot += fv; conn &= cv
        best_all = min(best_all, tot)
        if conn: best_conn = min(best_conn, tot)
    return best_all, best_conn

# ------------------------------------------------------------------ DP
def canon(blocks_of_bag):
    """blocks_of_bag: tuple of labels per bag vertex (bag sorted) -> canonical labels."""
    mp = {}; out = []
    for x in blocks_of_bag:
        if x not in mp: mp[x] = len(mp)
        out.append(mp[x])
    return tuple(out), mp

def dp_component(A, GL, comp, Z):
    sub = GL.subgraph(comp)
    if len(comp) == 1:
        return 0
    _, T = treewidth_min_fill_in(sub)
    root = next(iter(T.nodes))
    def gb(bagverts, part):  # g of bag parts, per block
        pass
    def rec(t, parent):
        bag = tuple(sorted(t))
        # table: dict (labels tuple, sizes tuple) -> value ; start as 'None' = no children
        tab = None
        for ch in T.neighbors(t):
            if ch == parent: continue
            ct = rec(ch, t)
            cb = tuple(sorted(ch))
            ct = transform(ct, cb, bag)
            tab = ct if tab is None else join(tab, ct, bag)
        if tab is None:
            # leaf: introduce all bag vertices from empty
            tab = transform({((), ()): 0}, (), bag)
        return tab
    def forget_vertex(tab, bag, v):
        i = bag.index(v); nb = bag[:i] + bag[i + 1:]; out = {}
        for (lab, sz), val in tab.items():
            nl = lab[:i] + lab[i + 1:]
            used = sorted(set(nl)); remap = {x: k for k, x in enumerate(used)}
            # block of v disappears if v alone in its block: closed
            newlab = []; mp = {}
            for x in nl:
                if x not in mp: mp[x] = len(mp)
                newlab.append(mp[x])
            nsz = [None] * len(mp)
            for x, k in mp.items(): nsz[k] = sz[x]
            key = (tuple(newlab), tuple(nsz))
            if key not in out or val < out[key]: out[key] = val
        return out, nb
    def introduce_vertex(tab, bag, v):
        nb = tuple(sorted(bag + (v,))); out = {}
        j = nb.index(v)
        for (lab, sz), val in tab.items():
            nblk = len(sz)
            # into existing block b
            for b in list(range(nblk)) + [None]:
                if b is None:
                    s = 0; k = 0; newsz = sz + (1,); nlabv = nblk
                else:
                    s = sz[b]
                    if s + 1 > Z: continue
                    k = sum(1 for u, x in zip(bag, lab) if x == b and GL.has_edge(u, v) and A[u][v])
                    newsz = sz[:b] + (sz[b] + 1,) + sz[b + 1:]; nlabv = b
                nl = lab[:j] + (nlabv,) + lab[j:]
                # careful: labels refer to blocks; newsz index by label
                key, mp = canon(nl)
                nsz = [None] * len(mp)
                for x, kk in mp.items(): nsz[kk] = newsz[x]
                key = (key, tuple(nsz))
                nv = val + s - 2 * k
                if key not in out or nv < out[key]: out[key] = nv
        return out, nb
    def transform(tab, cb, bag):
        cur = tab; curb = cb
        for v in cb:
            if v not in bag: cur, curb = forget_vertex(cur, curb, v)
        for v in bag:
            if v not in curb: cur, curb = introduce_vertex(cur, curb, v)
        assert curb == bag
        return cur
    def join(t1, t2, bag):
        out = {}
        by2 = {}
        for (lab, sz), val in t2.items(): by2.setdefault(lab, []).append((sz, val))
        for (lab, sz1), v1 in t1.items():
            for sz2, v2 in by2.get(lab, []):
                nb = len(sz1); tot = v1 + v2; nsz = []; ok = True
                for b in range(nb):
                    verts = [u for u, x in zip(bag, lab) if x == b]
                    bs = len(verts)
                    eb = sum(1 for u, w in itertools.combinations(verts, 2) if GL.has_edge(u, w) and A[u][w])
                    gbv = bs * (bs - 1) // 2 - 2 * eb
                    x1 = sz1[b] - bs; x2 = sz2[b] - bs
                    tot += -gbv + x1 * x2
                    s = sz1[b] + sz2[b] - bs
                    if s > Z: ok = False
                    nsz.append(s)
                if ok:
                    key = (lab, tuple(nsz))
                    if key not in out or tot < out[key]: out[key] = tot
        return out
    tab = rec(root, None)
    rb = tuple(sorted(root)); cur = tab
    for v in rb: cur, rb = forget_vertex(cur, rb, v) if False else forget_vertex(cur, rb, v)
    return min(cur.values())

def dp_m(A, st):
    GL = st['GL']; tot = 0
    for comp in nx.connected_components(GL):
        tot += dp_component(A, GL, sorted(comp), st['Z'])
    return tot

def main(seed, ninst, nr):
    rng = random.Random(seed)
    stt = dict(n=0, pairs=0, crit_all=0, crit_conn=0, crit_dp=0, sandwich=0, dpbad=0, pruned_pairs=0, big=0, big_nonempty=0, tw_drop=0)
    ex = []
    for it in range(ninst):
        n = rng.choice(nr)
        A, H = L.gen_instance(n, rng)
        if not H or not L.is_pair_maximal(A, H): continue
        opt = L.opt_subset(A)
        for l in range(0, 4):
            if l > 3: break
            st = setup(A, H, l)
            m_all, m_conn = brute_m(A, st)
            sumK = sum(f_true(A, K) for K in st['comps'] if any(v in st['BIG'] for v in K))
            m_dp = dp_m(A, st)
            base = int(A.sum() // 2) + sumK
            truth = opt <= len(H) + l
            stt['pairs'] += 1
            if st['BIG']: stt['big_nonempty'] += 1
            for name, m in (('all', m_all), ('conn', m_conn), ('dp', m_dp)):
                if (base + m <= len(H) + l) != truth:
                    stt['crit_' + name] += 1; ex.append((name, A.tolist(), H, l, opt, base, m))
            if not (m_all <= m_dp <= m_conn): stt['sandwich'] += 1; ex.append(('sand', A.tolist(), H, l, m_all, m_dp, m_conn))
            # Gamma' vs Gamma^l edge counts
            nS = sum(1 for i, j in L.prs(n) if st['S'][i, j] and len(st['Dl'][i, j]) > l and i not in st['BIG'] and j not in st['BIG'])
            stt['pruned_pairs'] += nS
        stt['n'] += 1
    return stt, ex

if __name__ == "__main__":
    seed = int(sys.argv[1]); ninst = int(sys.argv[2]); nr = [int(x) for x in sys.argv[3].split(',')]
    s, ex = main(seed, ninst, nr)
    print(s)
    for e in ex[:3]: print(e)

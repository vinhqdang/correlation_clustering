"""Wave 6 (ML) library, own code.  Graph = list of adjacency bitmasks adj[v].
P3 = (a,b,c) with centre b.  soft = pair inside a member of H, hard = any other pair.
OPT by exact subset DP (3^n), exact for n <= 13."""
import itertools, random
from functools import lru_cache

def key(a, b): return (a, b) if a < b else (b, a)

def pc(x): return bin(x).count("1")

def edges_of(adj):
    n = len(adj)
    return {(a, b) for a in range(n) for b in range(a + 1, n) if adj[a] >> b & 1}

def from_edges(n, E):
    adj = [0] * n
    for a, b in E:
        adj[a] |= 1 << b; adj[b] |= 1 << a
    return adj

def flip(adj, S):
    adj = list(adj)
    for a, b in S:
        adj[a] ^= 1 << b; adj[b] ^= 1 << a
    return adj

def p3s(adj):
    n = len(adj); out = []
    for b in range(n):
        nb = [x for x in range(n) if adj[b] >> x & 1]
        for a, c in itertools.combinations(nb, 2):
            if not (adj[a] >> c & 1): out.append((a, b, c))
    return out

def pairs3(t):
    a, b, c = t
    return {key(a, b), key(b, c), key(a, c)}

def is_cluster_graph(adj):
    n = len(adj)
    for b in range(n):
        nb = adj[b]
        x = nb
        while x:
            a = (x & -x).bit_length() - 1; x &= x - 1
            # N[a] must equal N[b]
            if (adj[a] | 1 << a) != (adj[b] | 1 << b): return False
    return True

def opt(adj, want_part=False):
    """min cluster editing cost via DP over subsets. returns cost (and one optimal partition)."""
    n = len(adj); N = 1 << n
    e_in = [0] * N  # edges inside mask
    for m in range(1, N):
        v = (m & -m).bit_length() - 1
        r = m & (m - 1)
        e_in[m] = e_in[r] + pc(adj[v] & r)
    tot = e_in[N - 1]
    f = [0] * N
    for m in range(1, N):
        s = pc(m); f[m] = s * (s - 1) // 2 - 2 * e_in[m]
    INF = 10 ** 9
    best = [INF] * N; best[0] = 0; arg = [0] * N
    for m in range(1, N):
        low = m & -m; rest = m ^ low
        sub = rest; bv = INF; ba = 0
        while True:
            c = sub | low
            v = f[c] + best[m ^ c]
            if v < bv: bv = v; ba = c
            if sub == 0: break
            sub = (sub - 1) & rest
        best[m] = bv; arg[m] = ba
    if not want_part: return tot + best[N - 1]
    part = []; m = N - 1
    while m:
        c = arg[m]; part.append(c); m ^= c
    return tot + best[N - 1], part

def all_opt_partitions(adj):
    """all optimal partitions (list of lists of masks) -- via recursion on DP table"""
    n = len(adj); N = 1 << n
    e_in = [0] * N
    for m in range(1, N):
        v = (m & -m).bit_length() - 1; r = m & (m - 1)
        e_in[m] = e_in[r] + pc(adj[v] & r)
    f = [0] * N
    for m in range(1, N):
        s = pc(m); f[m] = s * (s - 1) // 2 - 2 * e_in[m]
    INF = 10 ** 9
    best = [INF] * N; best[0] = 0
    for m in range(1, N):
        low = m & -m; rest = m ^ low; sub = rest; bv = INF
        while True:
            c = sub | low; v = f[c] + best[m ^ c]
            if v < bv: bv = v
            if sub == 0: break
            sub = (sub - 1) & rest
        best[m] = bv
    out = []
    def rec(m, acc):
        if m == 0: out.append(list(acc)); return
        low = m & -m; rest = m ^ low; sub = rest
        while True:
            c = sub | low
            if f[c] + best[m ^ c] == best[m]:
                acc.append(c); rec(m ^ c, acc); acc.pop()
            if sub == 0: break
            sub = (sub - 1) & rest
    rec(N - 1, [])
    return out, e_in[N - 1] + best[N - 1]

def maximal_pair_packing(adj, rng):
    ts = p3s(adj); rng.shuffle(ts); H = []; used = set()
    for t in ts:
        ps = pairs3(t)
        if not (ps & used): H.append(t); used |= ps
    return H

def is_pair_maximal(adj, H):
    used = set()
    for t in H: used |= pairs3(t)
    return all(pairs3(t) & used for t in p3s(adj))

def soft_pairs(H):
    s = set()
    for t in H: s |= pairs3(t)
    return s

def I_components(H):
    """connected components of I(H) (members adjacent iff they share a vertex); list of lists of member indices"""
    m = len(H); seen = [False] * m; comps = []
    for i in range(m):
        if seen[i]: continue
        st = [i]; seen[i] = True; comp = []
        while st:
            x = st.pop(); comp.append(x)
            for y in range(m):
                if not seen[y] and set(H[x]) & set(H[y]): seen[y] = True; st.append(y)
        comps.append(sorted(comp))
    return comps

def Iedges(H):
    return {(i, j) for i in range(len(H)) for j in range(i + 1, len(H)) if set(H[i]) & set(H[j])}

def hard_components(adj, H):
    n = len(adj); soft = soft_pairs(H)
    comp = [-1] * n; k = 0
    for s in range(n):
        if comp[s] >= 0: continue
        comp[s] = k; st = [s]
        while st:
            v = st.pop()
            for u in range(n):
                if adj[v] >> u & 1 and key(u, v) not in soft and comp[u] < 0:
                    comp[u] = k; st.append(u)
        k += 1
    return comp, k

def hard_deg(adj, H):
    n = len(adj); soft = soft_pairs(H)
    return [sum(1 for u in range(n) if adj[v] >> u & 1 and key(u, v) not in soft) for v in range(n)]

def random_instance(n, rng, mode=None):
    """planted cluster graph + k random flips; returns adj"""
    mode = mode or rng.choice(["planted", "gnp"])
    if mode == "gnp":
        p = rng.choice([0.25, 0.4, 0.55])
        E = {(a, b) for a in range(n) for b in range(a + 1, n) if rng.random() < p}
    else:
        r = rng.randint(2, max(2, n // 2)); lab = [rng.randrange(r) for _ in range(n)]
        E = {(a, b) for a in range(n) for b in range(a + 1, n) if lab[a] == lab[b]}
        for _ in range(rng.randint(1, max(2, n // 2))):
            a, b = rng.sample(range(n), 2); E ^= {key(a, b)}
    return from_edges(n, E)

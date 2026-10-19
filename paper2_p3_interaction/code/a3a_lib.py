"""Wave-3 audit library, written from scratch (does not import any earlier script)."""
import itertools, random
import numpy as np

def pidx(n):
    d = {}
    for k, (i, j) in enumerate(itertools.combinations(range(n), 2)):
        d[(i, j)] = d[(j, i)] = k
    return d

_rgs_cache = {}
def rgs(n):
    """all set partitions of range(n) as label arrays (restricted growth strings)."""
    if n in _rgs_cache: return _rgs_cache[n]
    out = []
    def rec(pref, m):
        if len(pref) == n:
            out.append(tuple(pref)); return
        for c in range(m + 1):
            rec(pref + [c], max(m, c + 1))
    rec([0], 1) if n > 0 else None
    arr = np.array(out, dtype=np.int8)
    _rgs_cache[n] = arr
    return arr

_same_cache = {}
def same_matrix(n):
    if n in _same_cache: return _same_cache[n]
    L = rgs(n)
    pr = list(itertools.combinations(range(n), 2))
    S = np.stack([(L[:, i] == L[:, j]) for i, j in pr], axis=1)
    _same_cache[n] = S
    return S

def costs_all(A):
    n = len(A)
    S = same_matrix(n).astype(np.int32)
    pr = list(itertools.combinations(range(n), 2))
    ed = np.array([A[i, j] for i, j in pr], dtype=np.int32)
    return S @ (1 - ed) + (1 - S) @ ed

def cost_labels(A, lab):
    n = len(A); c = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (lab[i] == lab[j]) != bool(A[i, j]): c += 1
    return c

def p3_triples(A):
    n = len(A); out = []
    for b in range(n):
        nb = [x for x in range(n) if A[b, x]]
        for a, c in itertools.combinations(nb, 2):
            if not A[a, c]: out.append((a, b, c))
    return out

def tri_pairs(t):
    a, b, c = t
    return [frozenset((a, b)), frozenset((b, c)), frozenset((a, c))]

def greedy_packing(A, rng, key=None, start=()):
    """random greedy PAIR-disjoint packing of induced P3s, maximal w.r.t. pair-disjointness."""
    T = p3_triples(A)
    rng.shuffle(T)
    if key: T.sort(key=key)
    used = set(); H = []
    for t in start:
        H.append(t); used.update(tri_pairs(t))
    for t in T:
        ps = tri_pairs(t)
        if all(p not in used for p in ps):
            H.append(t); used.update(ps)
    return H

def vertex_greedy_packing(A, rng):
    """maximal VERTEX-disjoint packing (weaker maximality)."""
    T = p3_triples(A); rng.shuffle(T)
    used = set(); H = []
    for t in T:
        if not (set(t) & used):
            H.append(t); used |= set(t)
    return H

def is_pair_maximal(A, H):
    used = set()
    for t in H: used.update(tri_pairs(t))
    return all(any(p in used for p in tri_pairs(t)) for t in p3_triples(A))

def info(A, H):
    n = len(A)
    VH = set(v for t in H for v in t)
    R = [v for v in range(n) if v not in VH]
    Q = {v: set() for v in range(n)}
    for t in H:
        for x in t: Q[x] |= set(t) - {x}
    # classes = components of G[R]
    seen = set(); classes = []
    for r in R:
        if r in seen: continue
        comp = [r]; seen.add(r); st = [r]
        while st:
            u = st.pop()
            for v in R:
                if v not in seen and A[u, v]: seen.add(v); comp.append(v); st.append(v)
        classes.append(sorted(comp))
    return VH, R, Q, classes

def lemma4_holds(A, R):
    for r, s in itertools.combinations(R, 2):
        n = len(A)
        if A[r, s]:
            for z in range(n):
                if z not in (r, s) and A[r, z] != A[s, z]: return False
        else:
            for z in range(n):
                if A[r, z] and A[s, z]: return False
    return True

def true_twin_pairs(A):
    n = len(A); out = []
    for u, v in itertools.combinations(range(n), 2):
        if A[u, v] and all(A[u, z] == A[v, z] for z in range(n) if z not in (u, v)): out.append((u, v))
    return out

# ---------------- generators ----------------
def sym(n): return np.zeros((n, n), dtype=np.int8)
def add(A, i, j, v=1): A[i, j] = A[j, i] = v

def gen_planted(rng, n, mode):
    k = rng.randint(2, max(2, n - 2))
    A = sym(n)
    p = rng.choice([.3, .5, .7])
    for i, j in itertools.combinations(range(k), 2):
        if rng.random() < p: add(A, i, j)
    v = k; classes = []
    while v < n:
        s = min(rng.randint(1, 4), n - v)
        cl = list(range(v, v + s)); v += s; classes.append(cl)
        for a, b in itertools.combinations(cl, 2): add(A, a, b)
    pool = list(range(k))
    used = set()
    hub = rng.randrange(k)
    for cl in classes:
        m = rng.randint(1, min(3, k))
        if mode == 'disjoint':
            avail = [z for z in pool if z not in used]
            att = rng.sample(avail, min(m, len(avail))) if avail else []
            used |= set(att)
        elif mode == 'hub':
            att = [hub] + ([rng.randrange(k)] if rng.random() < .3 else [])
        else:
            att = rng.sample(pool, m)
        for z in att:
            for x in cl: add(A, x, z)
    if rng.random() < .3:
        for _ in range(rng.randint(1, 2)):
            i, j = rng.sample(range(n), 2); add(A, i, j, 1 - A[i, j])
    return A, k

def gen_blowup(rng, n):
    m = rng.randint(3, max(3, n // 2))
    B = sym(m)
    for i, j in itertools.combinations(range(m), 2):
        if rng.random() < .5: add(B, i, j)
    sizes = [1] * m
    for _ in range(n - m): sizes[rng.randrange(m)] += 1
    lab = [i for i in range(m) for _ in range(sizes[i])]
    cliq = [rng.random() < .7 for _ in range(m)]
    A = sym(n)
    for u, v in itertools.combinations(range(n), 2):
        if lab[u] == lab[v]: A[u, v] = A[v, u] = int(cliq[lab[u]])
        else: A[u, v] = A[v, u] = B[lab[u], lab[v]]
    return A

def gen_gnp(rng, n):
    p = rng.choice([.25, .4, .55]); A = sym(n)
    for i, j in itertools.combinations(range(n), 2):
        if rng.random() < p: add(A, i, j)
    return A

def make_instance(rng, n):
    f = rng.choice(['pl_dis', 'pl_hub', 'pl_any', 'blow', 'gnp', 'pl_dis'])
    k = None
    if f.startswith('pl'):
        A, k = gen_planted(rng, n, f[3:]) if False else gen_planted(rng, n, {'dis': 'disjoint', 'hub': 'hub', 'any': 'any'}[f[3:]])
    elif f == 'blow': A = gen_blowup(rng, n)
    else: A = gen_gnp(rng, n)
    # packing: prefer triples inside core (to keep R large) or random
    pref = rng.choice(['core', 'rand', 'rand'])
    if k is not None and pref == 'core':
        key = lambda t: sum(1 for x in t if x >= k)
        H = greedy_packing(A, rng, key)
    else:
        H = greedy_packing(A, rng)
    return f, A, H

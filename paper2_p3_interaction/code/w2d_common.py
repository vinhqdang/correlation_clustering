"""Exact CEaMP tools (wave 2 search). ILP over final adjacency y with transitivity cuts (HiGHS via scipy.milp).
Graph = (n, frozenset of (a,b) a<b). Packing H = list of (a,b,c) induced paths a-b-c."""
import itertools, random
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix

def pid(n):
    return {p: i for i, p in enumerate(itertools.combinations(range(n), 2))}
def key(a, b): return (a, b) if a < b else (b, a)

_TR = {}
def trans(n):
    if n not in _TR:
        P = pid(n); r = []; c = []; v = []; k = 0
        for a, b, d in itertools.combinations(range(n), 3):
            ab, bd, ad = P[(a, b)], P[(b, d)], P[(a, d)]
            for x, y, z in ((ab, bd, ad), (ab, ad, bd), (bd, ad, ab)):  # y_x+y_y-y_z<=1
                r += [k, k, k]; c += [x, y, z]; v += [1, 1, -1]; k += 1
        _TR[n] = coo_matrix((v, (r, c)), shape=(k, n * (n - 1) // 2)).tocsr()
    return _TR[n]

def adj_vec(n, E):
    return np.array([1 if p in E else 0 for p in itertools.combinations(range(n), 2)])

def solve(n, E, w=None, cap=None, extra=None):
    """min sum_p cost(p) where edit var e_p=|y_p-g_p|. w: optional per-pair extra weight (small). cap: max edits.
    extra: list of (row vector over e, lo, hi). Returns (edits count, edit set) or None."""
    g = adj_vec(n, E); m = len(g)
    # y = g + s*e, s = 1-2g  => e = |y-g|;  cost = sum e. Work directly in e: y=g+s*e
    s = 1 - 2 * g
    T = trans(n)  # T y <= 1 -> T diag(s) e <= 1 - T g
    A = T.multiply(s[None, :]).tocsr()
    cons = [LinearConstraint(A, -np.inf, 1 - T @ g)]
    obj = np.ones(m) * 1000.0
    if w is not None: obj = obj + w
    if cap is not None: cons.append(LinearConstraint(np.ones((1, m)), 0, cap))
    for row, lo, hi in (extra or []): cons.append(LinearConstraint(np.asarray(row, float).reshape(1, -1), lo, hi))
    r = milp(obj, constraints=cons, integrality=np.ones(m), bounds=Bounds(0, 1))
    if r.status != 0: return None
    e = np.round(r.x).astype(int)
    ps = list(itertools.combinations(range(n), 2))
    S = frozenset(ps[i] for i in range(m) if e[i])
    return len(S), S

def opt(n, E):
    return solve(n, E)[0]

def induced_p3s(n, E):
    out = []
    for b in range(n):
        for a, c in itertools.combinations(range(n), 2):
            if a != b and c != b and key(a, b) in E and key(b, c) in E and key(a, c) not in E:
                out.append((a, b, c))
    return out
def ppairs(P):
    a, b, c = P; return [key(a, b), key(b, c), key(a, c)]

def rand_maximal_packing(rng, n, E):
    cand = induced_p3s(n, E); rng.shuffle(cand); used = set(); H = []
    for P in cand:
        ps = ppairs(P)
        if not any(q in used for q in ps): H.append(P); used.update(ps)
    return H

def max_packing_size(n, E):
    """exact maximum modification-disjoint induced-P3 packing (set packing ILP)."""
    cand = induced_p3s(n, E)
    if not cand: return 0, []
    P = pid(n); rows = []; cols = []
    for j, Q in enumerate(cand):
        for q in ppairs(Q): rows.append(P[q]); cols.append(j)
    A = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(P), len(cand))).tocsr()
    r = milp(-np.ones(len(cand)), constraints=LinearConstraint(A, 0, 1), integrality=np.ones(len(cand)), bounds=Bounds(0, 1))
    x = np.round(r.x).astype(int)
    return int(x.sum()), [cand[j] for j in range(len(cand)) if x[j]]

def vh(H): return {v for P in H for v in P}
def packed_pairs(H): return {q for P in H for q in ppairs(P)}

# ---------- graph families ----------
def gnp(rng, n, p): return frozenset(q for q in itertools.combinations(range(n), 2) if rng.random() < p)
def planted(rng, n, p_flip):
    k = rng.randint(2, max(2, n // 2)); lab = [rng.randrange(k) for _ in range(n)]
    return frozenset(q for q in itertools.combinations(range(n), 2) if (lab[q[0]] == lab[q[1]]) != (rng.random() < p_flip))
def blowup(rng, n, p):
    k = rng.randint(3, max(3, n // 2)); base = gnp(rng, k, p); sz = [1] * k
    for _ in range(n - k): sz[rng.randrange(k)] += 1
    cls = [i for i in range(k) for _ in range(sz[i])]; clique = [rng.random() < .6 for _ in range(k)]
    return frozenset(q for q in itertools.combinations(range(n), 2)
                     if (cls[q[0]] == cls[q[1]] and clique[cls[q[0]]]) or (cls[q[0]] != cls[q[1]] and key(cls[q[0]], cls[q[1]]) in base))
def near_clique(rng, n, p):
    k = rng.randint(n // 2, n); E = {q for q in itertools.combinations(range(k), 2) if rng.random() > p}
    E |= {(a, b) for a in range(k) for b in range(k, n) if rng.random() < .2}
    return frozenset(E)
def bip(rng, n, p):
    k = rng.randint(2, n - 2)
    return frozenset((a, b) for a in range(k) for b in range(k, n) if rng.random() < p)
def path_cliques(rng, n, p):
    sizes = []; r = n
    while r > 0: s = min(r, rng.randint(1, 4)); sizes.append(s); r -= s
    blocks = []; i = 0
    for s in sizes: blocks.append(list(range(i, i + s))); i += s
    E = set()
    for B in blocks: E |= set(itertools.combinations(B, 2))
    for A, B in zip(blocks, blocks[1:]):
        for a in A:
            for b in B:
                if rng.random() < p: E.add((a, b))
    if len(blocks) > 2 and rng.random() < .3:
        for a in blocks[0]:
            for b in blocks[-1]:
                if rng.random() < p: E.add((a, b))
    return frozenset(E)
FAMILIES = {"gnp": lambda r, n: gnp(r, n, r.choice([.2, .35, .5, .65])), "planted": lambda r, n: planted(r, n, r.choice([.05, .1, .2])),
            "blowup": lambda r, n: blowup(r, n, r.choice([.3, .5])), "nearclique": lambda r, n: near_clique(r, n, r.choice([.1, .2, .3])),
            "bip": lambda r, n: bip(r, n, r.choice([.4, .7, .9])), "pathcliques": lambda r, n: path_cliques(r, n, r.choice([.15, .4, .8]))}

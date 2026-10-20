"""Audit-10 library, written independently of w7_lib (numpy only).
Graph: n x n bool numpy adjacency. Member: ordered triple (a,b,c), edges ab, bc, non-edge ac. States 0=edit ab, 1=edit bc, 2=edit ac.
Ground truth for zero-excess (ZE) solutions: enumerate ALL set partitions of V (restricted growth strings), cost = #disagreeing pairs;
OPT>=|H| for a valid packing, ZE solutions = partitions of cost |H|; state vector read off from the edit set. No use of the decomposition lemma
except as an assertion (each ZE edit set must hit every member exactly once)."""
import itertools, os, random
import numpy as np

_PART = {}
def partitions(n):
    if n in _PART: return _PART[n]
    out = []
    def rec(pref, m):
        if len(pref) == n:
            out.append(pref[:]); return
        for l in range(m + 1):
            pref.append(l); rec(pref, max(m, l + 1)); pref.pop()
    rec([], 0)
    _PART[n] = np.array(out, dtype=np.int8)
    return _PART[n]

def pairs_of(t):
    a, b, c = t
    return [frozenset((a, b)), frozenset((b, c)), frozenset((a, c))]

def valid(A, H):
    seen = set()
    for (a, b, c) in H:
        if len({a, b, c}) < 3 or not (A[a, b] and A[b, c] and not A[a, c]): return False
        for p in pairs_of((a, b, c)):
            if p in seen: return False
            seen.add(p)
    return True

def owner_map(H):
    o = {}
    for i, t in enumerate(H):
        for p in pairs_of(t): o[p] = i
    return o

def soft_triangles(n, H):
    o = owner_map(H); res = []
    for p, q, r in itertools.combinations(range(n), 3):
        ow = [o.get(frozenset((p, q))), o.get(frozenset((q, r))), o.get(frozenset((p, r)))]
        if None not in ow and len(set(ow)) == 3: res.append((p, q, r))
    return res

def cmax(n, H):
    cnt = [0] * n
    for t in H:
        for v in t: cnt[v] += 1
    return max(cnt) if cnt else 0

def sigma_bruteforce(A, H):
    """set of ZE state vectors (tuples in {0,1,2}^m) from full partition enumeration"""
    n = len(A); m = len(H)
    P = partitions(n)
    same = P[:, :, None] == P[:, None, :]
    diff = same != A[None]
    iu = np.triu_indices(n, 1)
    D = diff[:, iu[0], iu[1]]           # N x npairs
    cost = D.sum(1)
    assert cost.min() >= m, "packing bound violated?"
    idx = {(int(i), int(j)): k for k, (i, j) in enumerate(zip(*iu))}
    def pk(x, y): return idx[(min(x, y), max(x, y))]
    mem = [[pk(a, b), pk(b, c), pk(a, c)] for (a, b, c) in H]
    sig = set()
    for r in np.nonzero(cost == m)[0]:
        st = []
        for ms in mem:
            e = [j for j in range(3) if D[r, ms[j]]]
            assert len(e) == 1
            st.append(e[0])
        sig.add(tuple(st))
    # sanity: total edits == m means no hard edit
    return sig, int(cost.min())

def opt_cost(A):
    n = len(A); P = partitions(n)
    same = P[:, :, None] == P[:, None, :]
    iu = np.triu_indices(n, 1)
    return int((same != A[None])[:, iu[0], iu[1]].sum(1).min())

def closed_all_rules(S, cap=40, rng=None, samples=20000):
    """closure of a set of state tuples under ALL 729 uniform majority rules (explicit enumeration of rule assignments per all-distinct pattern)"""
    S = list(S); Sset = set(S)
    if len(S) < 3: return True
    if len(S) > cap:
        rng = rng or random.Random(0)
        trip = [(rng.choice(S), rng.choice(S), rng.choice(S)) for _ in range(samples)]
    else:
        trip = itertools.product(S, S, S)
    for u, v, w in trip:
        fixed = []; var = {}
        for k, (x, y, z) in enumerate(zip(u, v, w)):
            if x == y or x == z: fixed.append(x)
            elif y == z: fixed.append(y)
            else:
                fixed.append(None); var.setdefault((x, y, z), []).append(k)
        if not var:
            if tuple(fixed) not in Sset: return False
            continue
        keys = list(var)
        for ch in itertools.product(range(3), repeat=len(keys)):
            res = fixed[:]
            for kk, c in zip(keys, ch):
                for pos in var[kk]: res[pos] = c
            if tuple(res) not in Sset: return False
    return True

def two_decomposable(S, m):
    S = set(S)
    if not S: return True
    proj = {(i, j): {(v[i], v[j]) for v in S} for i in range(m) for j in range(i + 1, m)}
    cnt = 0
    def rec(pref):
        nonlocal cnt
        i = len(pref)
        if i == m:
            cnt += 1; return
        for s in range(3):
            if all((pref[j], s) in proj[(j, i)] for j in range(i)): rec(pref + [s])
    rec([])
    return cnt == len(S)

def planted(n, m, rng, big=True):
    """random planted tight-ish instance: cluster graph C, members = triples with one flipped pair making an induced P3. Returns (A,H) or None"""
    k = rng.randint(2, max(2, n // 2)) if not big else rng.randint(2, 4)
    lab = [rng.randrange(k) for _ in range(n)]
    A = np.zeros((n, n), dtype=bool)
    for i in range(n):
        for j in range(i + 1, n):
            if lab[i] == lab[j]: A[i, j] = A[j, i] = True
    used = set(); H = []; flips = []
    for _ in range(60):
        if len(H) >= m: break
        a, b, c = rng.sample(range(n), 3)
        same = [lab[a] == lab[b], lab[b] == lab[c], lab[a] == lab[c]]
        s = sum(same)
        if s == 3:
            # flip one pair, centre = the vertex opposite
            f = rng.choice([(a, b), (b, c), (a, c)]); x, y = f
            z = ({a, b, c} - {x, y}).pop()
            tri = (x, z, y)       # centre z
        elif s == 1:
            # two in the same cluster, third apart: flip a cross pair to make it an edge
            if same[0]: x, y, o = a, b, c
            elif same[1]: x, y, o = b, c, a
            else: x, y, o = a, c, b
            t = rng.choice([x, y]); f = (t, o)
            other = y if t == x else x
            tri = (other, t, o)   # centre t
        else: continue
        ps = [frozenset(p) for p in itertools.combinations((a, b, c), 2)]
        if any(p in used for p in ps): continue
        used.update(ps); flips.append(f); H.append(tri)
    if len(H) < m: return None
    for x, y in flips: A[x, y] = A[y, x] = not A[x, y]
    if not valid(A, H): return None
    return A, H

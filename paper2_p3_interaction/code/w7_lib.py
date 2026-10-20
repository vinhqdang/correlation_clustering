"""Wave 7 library (written from scratch; numpy + networkx only).
Vertices 0..n-1; graph = symmetric bool matrix A; member = (a,b,c) with centre b, edges ab, bc, non-edge ac.
Zero-excess (ZE) solutions: by the decomposition lemma a solution of cost |H| edits exactly one pair of every
member and nothing else, so Sigma = {state vectors in {L,R,M}^H : G xor edits is a cluster graph}; L=edit ab, R=edit bc, M=edit ac.
Tightness (OPT = |H|) <=> Sigma non-empty (OPT >= |H| always for a valid packing)."""
import itertools, random
import numpy as np
import networkx as nx

ST = 'LRM'

def key(a, b):
    return (a, b) if a < b else (b, a)

def mpairs(t):
    a, b, c = t
    return [key(a, b), key(b, c), key(a, c)]

class Inst:
    def __init__(self, A, H):
        self.A = np.array(A, dtype=bool)
        self.n = len(self.A)
        self.H = [tuple(t) for t in H]
        self.owner = {}
        for i, t in enumerate(self.H):
            for p in mpairs(t):
                self.owner[p] = i

    def valid(self):
        A = self.A
        seen = set()
        for (a, b, c) in self.H:
            if len({a, b, c}) < 3:
                return False
            if not (A[a, b] and A[b, c] and not A[a, c]):
                return False
            for p in mpairs((a, b, c)):
                if p in seen:
                    return False
                seen.add(p)
        return True

    def cmax(self):
        cnt = [0] * self.n
        for t in self.H:
            for v in t:
                cnt[v] += 1
        return max(cnt) if cnt else 0

    def Ig(self):
        g = nx.Graph()
        g.add_nodes_from(range(len(self.H)))
        for i in range(len(self.H)):
            for j in range(i + 1, len(self.H)):
                if set(self.H[i]) & set(self.H[j]):
                    g.add_edge(i, j)
        return g

    def soft_triangles(self):
        """vertex triples whose three pairs are packed by three distinct members"""
        res = []
        for p, q, r in itertools.combinations(range(self.n), 3):
            o = [self.owner.get(key(p, q)), self.owner.get(key(q, r)), self.owner.get(key(p, r))]
            if None not in o and len(set(o)) == 3:
                res.append((p, q, r))
        return res

    def hard_deg(self, v):
        return sum(1 for x in range(self.n) if x != v and self.A[v, x] and key(v, x) not in self.owner)

    def sigma(self):
        """all ZE state vectors as list of strings"""
        m = len(self.H)
        n = self.n
        if m == 0:
            return ['']
        vecs = np.array(list(itertools.product(range(3), repeat=m)), dtype=np.int8)  # (V,m)
        V = len(vecs)
        ed = np.zeros((m, 3, 2), dtype=np.int64)
        for i, (a, b, c) in enumerate(self.H):
            ed[i, 0] = (a, b)
            ed[i, 1] = (b, c)
            ed[i, 2] = (a, c)
        Ab = np.broadcast_to(self.A, (V, n, n)).copy()
        for i in range(m):
            u = ed[i, vecs[:, i], 0]
            v = ed[i, vecs[:, i], 1]
            idx = np.arange(V)
            Ab[idx, u, v] ^= True
            Ab[idx, v, u] ^= True
        B = Ab | np.eye(n, dtype=bool)[None]
        Bf = B.astype(np.float32)
        P = np.matmul(Bf, Bf) > 0.5
        ok = np.all(~P | B, axis=(1, 2))
        return [''.join(ST[x] for x in vecs[k]) for k in np.nonzero(ok)[0]]

def mmaj(a, b, c):
    """M-majority on strings of states"""
    out = []
    for x, y, z in zip(a, b, c):
        if x == y or x == z:
            out.append(x)
        elif y == z:
            out.append(y)
        else:
            out.append('M')
    return ''.join(out)

def closed_mmaj(vecs):
    vs = set(vecs)
    for a, b, c in itertools.combinations_with_replacement(vecs, 3):
        for t in {(a, b, c), (a, c, b), (b, c, a)}:
            pass
        if mmaj(a, b, c) not in vs:
            return False
    return True

def closed_mmaj_fast(vecs):
    vs = set(vecs)
    L = list(vecs)
    for a in L:
        for b in L:
            for c in L:
                if mmaj(a, b, c) not in vs:
                    return False
    return True

def two_decomposable(vecs, m=None):
    """True iff the set equals the set of vectors all of whose coordinate pairs lie in the binary projections."""
    if not vecs:
        return True
    m = len(vecs[0])
    vset = set(vecs)
    proj = {(i, j): {(v[i], v[j]) for v in vecs} for i in range(m) for j in range(i + 1, m)}
    sol = set()
    def rec(pref):
        i = len(pref)
        if i == m:
            sol.add(''.join(pref))
            return
        for s in ST:
            if all((pref[j], s) in proj[(j, i)] for j in range(i)):
                rec(pref + [s])
    rec([])
    return sol == vset

def gen_planted(n, m, rng, cmax=2, maxcl=5, tries=400, want_triangle=False, forest=False):
    """random planted tight instance. returns Inst or None."""
    lab = [rng.randrange(max(2, n // 2)) for _ in range(n)]
    cl = {}
    for v, l in enumerate(lab):
        cl.setdefault(l, []).append(v)
    if max(len(x) for x in cl.values()) > maxcl:
        return None
    A = np.zeros((n, n), dtype=bool)
    for l, vs in cl.items():
        for u, v in itertools.combinations(vs, 2):
            A[u, v] = A[v, u] = True
    used = set()
    cnt = [0] * n
    H = []
    edits = []
    for _ in range(tries):
        if len(H) >= m:
            break
        s = rng.choice(ST)
        if s == 'L':       # a apart from b,c ; b,c together
            b, c = rng.sample(range(n), 2)
            if lab[b] != lab[c]:
                continue
            a = rng.randrange(n)
            if lab[a] == lab[b]:
                continue
        elif s == 'R':
            a, b = rng.sample(range(n), 2)
            if lab[a] != lab[b]:
                continue
            c = rng.randrange(n)
            if lab[c] == lab[a]:
                continue
        else:
            a, b, c = rng.sample(range(n), 3)
            if not (lab[a] == lab[b] == lab[c]):
                continue
        if len({a, b, c}) < 3:
            continue
        ps = mpairs((a, b, c))
        if any(p in used for p in ps):
            continue
        if any(cnt[v] >= cmax for v in (a, b, c)):
            continue
        # set the soft pairs to P3 pattern
        if s == 'L':
            edits.append(key(a, b))
        elif s == 'R':
            edits.append(key(b, c))
        else:
            edits.append(key(a, c))
        for p in ps:
            used.add(p)
        for v in (a, b, c):
            cnt[v] += 1
        H.append((a, b, c))
    if len(H) < m:
        return None
    for (u, v) in edits:
        A[u, v] = A[v, u] = not A[u, v]
    I = Inst(A, H)
    if not I.valid():
        return None
    return I

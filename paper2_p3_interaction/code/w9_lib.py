"""Wave 9 library (written fresh; numpy/itertools only).
Vertices 0..n-1, E = set of edges (u<v).  Member (a,b,c): centre b, ab,bc in E, ac not in E.
Member pairs in order (ab, bc, ac); an 'edit mask' is a 3-bit subset (bit0=ab, bit1=bc, bit2=ac) of edited pairs.
Final adjacency of a soft pair after its member's edit: y_e = [e in E] xor [e edited].
A member's own triple is a valid cluster graph for exactly 5 masks: L=001(edit ab), R=010(edit bc), M=100(edit ac),
O=011 (edit ab,bc: all apart), X=111 (edit all: a,c together, b apart).
cost(L)=cost(R)=cost(M)=1, cost(O)=2, cost(X)=3.
"""
import itertools, random, sys
import numpy as np
import networkx as nx

def key(a, b):
    return (a, b) if a < b else (b, a)

NAMES = {1: 'L', 2: 'R', 4: 'M', 3: 'O', 7: 'X'}
VALID = [1, 2, 4, 3, 7]
COST = {1: 1, 2: 1, 4: 1, 3: 2, 7: 3}
ZE = [1, 2, 4]

def mpairs(t):
    a, b, c = t
    return [key(a, b), key(b, c), key(a, c)]

class Inst:
    def __init__(self, n, E, H):
        self.n = n
        self.E = set(key(*e) for e in E)
        self.H = [tuple(t) for t in H]
        self.owner = {}
        for i, t in enumerate(self.H):
            for k, p in enumerate(mpairs(t)):
                assert p not in self.owner
                self.owner[p] = (i, k)
        self.adj = [[False] * n for _ in range(n)]
        for u, v in self.E:
            self.adj[u][v] = self.adj[v][u] = True
        self.triples = list(itertools.combinations(range(n), 3))

    def isE(self, u, v):
        return self.adj[u][v]

    def valid_packing(self):
        for (a, b, c) in self.H:
            if not (self.isE(a, b) and self.isE(b, c) and not self.isE(a, c)):
                return False
        return True

    def pair_maximal(self):
        for (p, q, r) in self.triples:
            e = [self.isE(p, q), self.isE(q, r), self.isE(p, r)]
            if sum(e) == 2 and not any(key(*x) in self.owner for x in ((p, q), (q, r), (p, r))):
                return False
        return True

    def Ig(self):
        g = nx.Graph()
        g.add_nodes_from(range(len(self.H)))
        for i in range(len(self.H)):
            for j in range(i + 1, len(self.H)):
                if set(self.H[i]) & set(self.H[j]):
                    g.add_edge(i, j)
        return g

    def cmax(self):
        cnt = [0] * self.n
        for t in self.H:
            for v in t:
                cnt[v] += 1
        return max(cnt) if cnt else 0

    def soft_triangles(self):
        res = []
        for (p, q, r) in self.triples:
            o = [self.owner.get(key(p, q)), self.owner.get(key(q, r)), self.owner.get(key(p, r))]
            if None not in o and len({x[0] for x in o}) == 3:
                res.append((p, q, r))
        return res

    # ---- exact CSP evaluation -------------------------------------------------
    def compile(self, U=()):
        """triple constraint list for G xor U: each entry (owners tuple, fn) ; owners = distinct member indices with a soft pair in the triple;
        constraint evaluated from the masks of those owners."""
        Uset = set(key(*u) for u in U)
        cons = []
        for (p, q, r) in self.triples:
            prs = [key(p, q), key(q, r), key(p, r)]
            base = []
            soft = []
            for k, e in enumerate(prs):
                a = self.isE(*e) ^ (e in Uset)
                base.append(a)
                if e in self.owner:
                    soft.append((k, self.owner[e]))
            if not soft:
                if sum(base) == 2:
                    cons.append(((), None, base, soft))  # static violation
                continue
            owners = tuple(sorted({o[0] for _, o in soft}))
            cons.append((owners, None, base, soft))
        return cons

    def feasible_triple(self, base, soft, masks):
        cnt = 0
        for k in range(3):
            v = base[k]
            for kk, (i, pk) in soft:
                pass
        return None

def triple_ok(base, soft, masks):
    """base: three base adjacencies (after U); soft: list of (k, (member, pairindex)); masks: dict member->mask"""
    s = 0
    sd = {k: o for k, o in soft}
    for k in range(3):
        a = base[k]
        if k in sd:
            i, pk = sd[k]
            if masks[i] >> pk & 1:
                a = not a
        s += a
    return s != 2

class Solver:
    """exact min cost of soft edits (members' masks) given hard edits U (set of hard pairs flipped), by DFS over members"""
    def __init__(self, inst, U=(), order=None, allowed=None):
        self.inst = inst
        self.U = set(key(*u) for u in U)
        m = len(inst.H)
        self.m = m
        self.order = list(range(m)) if order is None else order
        self.allowed = allowed  # dict member -> list of masks allowed (default VALID)
        cons = inst.compile(self.U)
        self.static_bad = any(len(c[0]) == 0 for c in cons)
        pos = {mm: k for k, mm in enumerate(self.order)}
        self.by_last = [[] for _ in range(m)]
        for owners, _, base, soft in cons:
            if not owners:
                continue
            last = max(pos[o] for o in owners)
            self.by_last[last].append((owners, base, soft))
        # reduce: constraints with only own-member triple pairs are handled by VALID; keep anyway
    def solve(self, bound=10 ** 9, want_all=False):
        if self.static_bad:
            return None
        m = self.m
        best = [bound + 1]
        sols = []
        masks = {}
        order = self.order
        allowed = self.allowed
        def rec(k, cost):
            if cost >= best[0] and not (want_all and cost == best[0]):
                return
            if k == m:
                if cost < best[0]:
                    best[0] = cost
                    sols.clear()
                if want_all or not sols:
                    sols.append(dict(masks))
                return
            i = order[k]
            opts = VALID if allowed is None or i not in allowed else allowed[i]
            for mk in sorted(opts, key=lambda x: COST[x]):
                masks[i] = mk
                ok = True
                for owners, base, soft in self.by_last[k]:
                    if not triple_ok(base, soft, masks):
                        ok = False
                        break
                if ok:
                    rec(k + 1, cost + COST[mk])
            masks.pop(i, None)
        rec(0, 0)
        if best[0] > bound:
            return None
        return best[0], sols

def m_of(inst, U=(), bound=10 ** 9):
    r = Solver(inst, U).solve(bound)
    return None if r is None else r[0]

# ---------------------------------------------------------------------------
def maximal_packing(n, E, rng, forest=True, cmax=2):
    """random pair-maximal modification-disjoint packing (greedy over shuffled induced P3s).  returns H or None"""
    E = set(E)
    adj = [[False] * n for _ in range(n)]
    for u, v in E:
        adj[u][v] = adj[v][u] = True
    p3 = []
    for b in range(n):
        nb = [x for x in range(n) if adj[b][x]]
        for x, y in itertools.combinations(nb, 2):
            if not adj[x][y]:
                p3.append((x, b, y))
    rng.shuffle(p3)
    used = set()
    H = []
    cnt = [0] * n
    for (a, b, c) in p3:
        ps = [key(a, b), key(b, c), key(a, c)]
        if any(p in used for p in ps):
            continue
        if any(cnt[v] >= cmax for v in (a, b, c)):
            continue
        H.append((a, b, c))
        for p in ps:
            used.add(p)
        for v in (a, b, c):
            cnt[v] += 1
    return H

def gen_instance(rng, nmin=8, nmax=13, cmax=2, forest=True, noise=None, tries=50, maxcl=5):
    for _ in range(tries):
        n = rng.randint(nmin, nmax)
        k = rng.randint(2, max(2, n // 2))
        lab = [rng.randrange(k) for _ in range(n)]
        E = set()
        for u, v in itertools.combinations(range(n), 2):
            same = lab[u] == lab[v]
            pr = noise if noise is not None else rng.choice([0.05, 0.1, 0.15, 0.25])
            if same:
                if rng.random() > pr:
                    E.add((u, v))
            else:
                if rng.random() < pr * 0.6:
                    E.add((u, v))
        H = maximal_packing(n, E, rng, cmax=cmax)
        if not H:
            continue
        inst = Inst(n, E, H)
        if not inst.pair_maximal():
            continue
        if forest and not nx.is_forest(inst.Ig()):
            continue
        if inst.soft_triangles():
            continue
        if len(H) > 8:
            continue
        return inst
    return None

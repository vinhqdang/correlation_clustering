"""w11_lib.py -- wave 11 library (written from scratch).

Instance: n vertices, edge set E, packing H = list of members (a,b,c): centre b, ab,bc in E, ac not in E.
Pair classes: soft (owned by a member), local hard (both ends covered, some member through u equals / shares a vertex with some member through v),
remote (everything else; soft pairs are local).  F'(T): min #edits over assignments to soft+local pairs, every member edited at least once,
such that every triple all of whose remote pairs lie in T is not an induced P3 (remote pairs keep their G-status; other triples dropped).
A budget t=|H|+l; T is a certificate iff F'(T) > t.

Oracles
  FPrime(inst, t)      incremental SAT (pysat) solver:  feasible(T) -> assignment or None   (F'(T) <= t ?)
  min_certificate(...)  exact minimum certificate by CEGAR (CP-SAT master over remote pairs, SAT checker)
  opt(inst)            exact OPT by subset DP (n <= 14) -- independent of everything else
"""
import itertools, random, time
from pysat.solvers import Cadical153, Glucose4
from pysat.card import ITotalizer
from ortools.sat.python import cp_model


def key(u, v):
    return (u, v) if u < v else (v, u)


class Inst:
    def __init__(self, n, E, H):
        self.n = n
        self.E = set(key(*e) for e in E)
        self.H = [tuple(h) for h in H]
        self.adj = [[False] * n for _ in range(n)]
        for u, v in self.E:
            self.adj[u][v] = self.adj[v][u] = True
        self.owner = {}
        for i, (a, b, c) in enumerate(self.H):
            for k, p in enumerate((key(a, b), key(b, c), key(a, c))):
                self.owner[p] = (i, k)
        self.mem = [[] for _ in range(n)]
        for i, t in enumerate(self.H):
            for v in t:
                self.mem[v].append(i)
        self._nm = {}
        self._cls = None

    # ---- validity -------------------------------------------------
    def valid(self):
        seen = set()
        for (a, b, c) in self.H:
            if not (self.adj[a][b] and self.adj[b][c] and not self.adj[a][c]):
                return False
            for p in (key(a, b), key(b, c), key(a, c)):
                if p in seen:
                    return False
                seen.add(p)
        return True

    def pair_maximal(self):
        for x, y, z in itertools.combinations(range(self.n), 3):
            ps = [key(x, y), key(y, z), key(x, z)]
            if any(p in self.owner for p in ps):
                continue
            if sum(self.adj[p[0]][p[1]] for p in ps) == 2:
                return False
        return True

    def cmax(self):
        return max((len(m) for m in self.mem), default=0)

    def nbrs(self, i):
        if i not in self._nm:
            s = set()
            for v in self.H[i]:
                s.update(self.mem[v])
            self._nm[i] = s
        return self._nm[i]

    def local(self, u, v):
        if not self.mem[u] or not self.mem[v]:
            return False
        for i in self.mem[u]:
            ni = self.nbrs(i)
            for j in self.mem[v]:
                if j in ni:
                    return True
        return False

    def classify(self):
        if self._cls is None:
            soft, lochard, remote = [], [], []
            for p in itertools.combinations(range(self.n), 2):
                if p in self.owner:
                    soft.append(p)
                elif self.local(*p):
                    lochard.append(p)
                else:
                    remote.append(p)
            self._cls = (soft, lochard, remote)
        return self._cls

    def remote(self):
        return self.classify()[2]

    def IH(self):
        import networkx as nx
        g = nx.Graph()
        g.add_nodes_from(range(len(self.H)))
        for i, j in itertools.combinations(range(len(self.H)), 2):
            if set(self.H[i]) & set(self.H[j]):
                g.add_edge(i, j)
        return g

    def soft_triangles(self):
        res = []
        for x, y, z in itertools.combinations(range(self.n), 3):
            o = [self.owner.get(key(x, y)), self.owner.get(key(y, z)), self.owner.get(key(x, z))]
            if None not in o and len({t[0] for t in o}) == 3:
                res.append((x, y, z))
        return res


# =====================================================================
class FPrime:
    """incremental oracle for F'(T) <= t.  Remote-pair selectors t_e are assumption literals."""

    def __init__(self, inst, t):
        self.I = inst
        self.t = t
        soft, lochard, remote = inst.classify()
        self.varpairs = soft + lochard
        self.remote = remote
        self.vid = {}
        nv = 0
        for p in self.varpairs:
            nv += 1
            self.vid[p] = nv           # e_p: pair edited
        self.sel = {}
        for p in remote:
            nv += 1
            self.sel[p] = nv           # t_e: remote pair in T
        clauses = []
        # member non-empty
        for (a, b, c) in inst.H:
            clauses.append([self.vid[key(a, b)], self.vid[key(b, c)], self.vid[key(a, c)]])
        # triples
        self.trip = []   # (pairs, remote-set) for every triple that contains a remote pair
        for x, y, z in itertools.combinations(range(inst.n), 3):
            ps = (key(x, y), key(y, z), key(x, z))
            rem = frozenset(p for p in ps if p in self.sel)
            vp = [p for p in ps if p in self.vid]
            self.trip.append((ps, rem))
            # forbid exactly-two-edges: for each bad pattern (one non-edge position) clause
            for zero in range(3):
                pat = [1, 1, 1]
                pat[zero] = 0
                cl = []
                sat = False
                guards = []
                for p, want in zip(ps, pat):
                    g = 1 if p in inst.E else 0
                    if p in self.vid:
                        # status s = g xor e ; clause literal true iff s != want, i.e. e = (g != (1-want)) ... compute:
                        # s != want  <=>  (g xor e) != want  <=> e != (g xor want)  <=> e == 1-(g xor want)
                        need_e = 1 - (g ^ want)
                        cl.append(self.vid[p] if need_e == 1 else -self.vid[p])
                    else:  # remote constant g
                        if g != want:
                            sat = True
                            break
                        guards.append(-self.sel[p])
                if sat:
                    continue
                clauses.append(cl + guards)
        self.clauses = clauses
        self.solver = Cadical153(bootstrap_with=clauses)
        lits = [self.vid[p] for p in self.varpairs]
        self.tot = ITotalizer(lits=lits, ubound=t + 1, top_id=nv)
        for cl in self.tot.cnf.clauses:
            self.solver.add_clause(cl)
        self.rhs = self.tot.rhs

    def feasible(self, T):
        """return dict pair->edited(0/1) for non-remote pairs if F'(T) <= t else None"""
        Ts = set(T)
        ass = [(self.sel[p] if p in Ts else -self.sel[p]) for p in self.remote]
        ass.append(-self.rhs[self.t])   # at most t true
        if self.solver.solve(assumptions=ass):
            m = self.solver.get_model()
            ms = set(l for l in m if l > 0)
            return {p: (1 if self.vid[p] in ms else 0) for p in self.varpairs}
        return None

    def violated_sets(self, edits):
        """remote-sets R(tau) of triples violated by the assignment (edits on var pairs, remote constant)."""
        I = self.I
        out = set()
        for ps, rem in self.trip:
            if not rem:
                continue
            s = 0
            for p in ps:
                g = 1 if p in I.E else 0
                if p in edits:
                    g ^= edits[p]
                s += g
            if s == 2:
                out.add(rem)
        return out

    def cost(self, edits):
        return sum(edits.values())


def minimalize(sets):
    sets = sorted(set(sets), key=len)
    keep = []
    for s in sets:
        if not any(k <= s for k in keep):
            keep.append(s)
    return keep


def min_certificate(inst, l, time_limit=600, verbose=False, init_T=None, max_size=None, workers=4):
    """exact minimum certificate for budget t=|H|+l. Returns (size, T) or (None,None) if F'(N) <= t (no certificate)."""
    t = len(inst.H) + l
    O = FPrime(inst, t)
    N = list(O.remote)
    if O.feasible(N) is not None:
        return None, None
    cons = []      # list of lists of remote-sets (violation families)
    cache = {}
    t0 = time.time()
    # seed: a few random models
    rnd = random.Random(1)
    best_lb = 0
    it = 0
    T = set(init_T or [])
    while True:
        it += 1
        x = O.feasible(T)
        if x is None:
            return len(T), sorted(T)
        fam = minimalize(O.violated_sets(x))
        cons.append(fam)
        # master: min sum t_e s.t. for every fam there is R in fam with all t_e (e in R) true
        m = cp_model.CpModel()
        tv = {p: m.NewBoolVar('t%d_%d' % p) for p in N}
        ys = {}
        for fam_ in cons:
            ors = []
            for R in fam_:
                if len(R) == 1:
                    ors.append(tv[next(iter(R))])
                else:
                    if R not in ys:
                        y = m.NewBoolVar('y%d' % len(ys))
                        for e in R:
                            m.AddImplication(y, tv[e])
                        ys[R] = y
                    ors.append(ys[R])
            m.AddBoolOr(ors)
        m.Minimize(sum(tv.values()))
        s = cp_model.CpSolver()
        s.parameters.num_workers = workers
        s.parameters.max_time_in_seconds = max(10, time_limit - (time.time() - t0))
        st = s.Solve(m)
        if st != cp_model.OPTIMAL:
            raise RuntimeError('master status %s' % s.StatusName(st))
        T = set(p for p in N if s.Value(tv[p]))
        if verbose and it % 10 == 0:
            print('  it', it, '|T|', len(T), 'cons', len(cons), round(time.time() - t0, 1), flush=True)
        if max_size is not None and len(T) > max_size:
            return len(T), None   # lower bound exceeded


def greedy_minimal(inst, l, order_seed=0):
    """inclusion-minimal certificate by deletion (upper bound on min certificate)."""
    t = len(inst.H) + l
    O = FPrime(inst, t)
    N = list(O.remote)
    if O.feasible(N) is not None:
        return None
    T = list(N)
    rnd = random.Random(order_seed)
    rnd.shuffle(T)
    cur = set(T)
    for e in T:
        cur.discard(e)
        if O.feasible(cur) is not None:
            cur.add(e)
    return sorted(cur)


# =====================================================================
def opt(inst):
    """exact OPT for n <= 16 by subset DP. cost = |E| + sum over clusters (non-edges inside - edges inside)."""
    n = inst.n
    assert n <= 16
    N = 1 << n
    w = [0] * N
    # w[S] = nonedges(S) - edges(S); compute incrementally
    for S in range(1, N):
        low = (S & -S).bit_length() - 1
        R = S & (S - 1)
        add = 0
        r = R
        while r:
            b = (r & -r).bit_length() - 1
            add += -1 if inst.adj[low][b] else 1
            r &= r - 1
        w[S] = w[R] + add
    INF = 10 ** 9
    f = [INF] * N
    f[0] = 0
    for S in range(1, N):
        low = S & -S
        rest = S ^ low
        best = INF
        sub = rest
        while True:
            C = sub | low
            v = w[C] + f[S ^ C]
            if v < best:
                best = v
            if sub == 0:
                break
            sub = (sub - 1) & rest
        f[S] = best
    return len(inst.E) + f[N - 1]


def fprime_value(inst, T, tmax=None):
    """exact F'(T) (min edits) by increasing budget; None if > tmax"""
    base = len(inst.H)
    t = base
    lim = tmax if tmax is not None else base + 12
    while t <= lim:
        O = FPrime(inst, t)
        if O.feasible(T) is not None:
            return t
        t += 1
    return None


def min_increment(inst, Phi, time_limit=600, max_size=None, workers=4, tabs=None):
    """min |C| (C subset of remote, disjoint from Phi) with F'(Phi u C) > F'(Phi). Returns (size, C, F'(Phi)) or (None,None,F'(Phi)) if F'(N)=F'(Phi)."""
    Phi = set(Phi)
    t = fprime_value(inst, Phi) if tabs is None else tabs
    O = FPrime(inst, t)
    N = list(O.remote)
    if O.feasible(N) is not None:
        return None, None, t
    cons = []
    T = set(Phi)
    t0 = time.time()
    while True:
        x = O.feasible(T)
        if x is None:
            return len(T - Phi), sorted(T - Phi), t
        cons.append(minimalize(O.violated_sets(x)))
        m = cp_model.CpModel()
        tv = {p: m.NewBoolVar('t%d_%d' % p) for p in N}
        for p in Phi: m.Add(tv[p] == 1)
        ys = {}
        for fam_ in cons:
            ors = []
            for R in fam_:
                if len(R) == 1:
                    ors.append(tv[next(iter(R))])
                else:
                    if R not in ys:
                        y = m.NewBoolVar('y%d' % len(ys))
                        for e in R: m.AddImplication(y, tv[e])
                        ys[R] = y
                    ors.append(ys[R])
            m.AddBoolOr(ors)
        m.Minimize(sum(tv[p] for p in N if p not in Phi))
        s = cp_model.CpSolver(); s.parameters.num_workers = workers
        s.parameters.max_time_in_seconds = max(10, time_limit - (time.time() - t0))
        st = s.Solve(m)
        if st != cp_model.OPTIMAL: raise RuntimeError('master %s' % s.StatusName(st))
        T = set(p for p in N if s.Value(tv[p]))
        if max_size is not None and len(T - Phi) > max_size:
            return len(T - Phi), None, t

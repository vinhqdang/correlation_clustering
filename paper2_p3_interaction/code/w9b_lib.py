"""w9b_lib.py -- wave 9b library, independent of earlier waves.
Instances (n, E, members); members are (a,b,c) with centre b, ab, bc in E, ac not in E.
Two independent oracles:
  * opt_cpsat(inst, ...)   exact minimum cluster-editing cost (CP-SAT over pair variables with transitivity)  [optionally with constraints]
  * ze_states(inst)        all zero-excess state vectors (L=0,R=1,M=2) by DFS over members with cluster-graph test on the full vertex set
Locality: pair {u,v} is LOCAL iff both ends lie in V(H) and some member through u equals or shares a vertex with some member through v.
"""
import itertools, random
from ortools.sat.python import cp_model

L_, R_, M_ = 0, 1, 2

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
                self.owner[p] = (i, k)       # k: 0 = ab, 1 = bc, 2 = ac
        self.mem = [[] for _ in range(n)]
        for i, t in enumerate(self.H):
            for v in t:
                self.mem[v].append(i)

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

    def nbr_members(self, i):
        s = set()
        for v in self.H[i]:
            s.update(self.mem[v])
        return s

    def local(self, u, v):
        if not self.mem[u] or not self.mem[v]:
            return False
        for i in self.mem[u]:
            ni = self.nbr_members(i)
            for j in self.mem[v]:
                if j in ni:
                    return True
        return False

    def hard_pairs(self):
        return [(u, v) for u, v in itertools.combinations(range(self.n), 2) if (u, v) not in self.owner]

    def far_pairs(self):
        return [p for p in self.hard_pairs() if not self.local(*p)]

    def local_hard_pairs(self):
        return [p for p in self.hard_pairs() if self.local(*p)]

    def IH(self):
        import networkx as nx
        g = nx.Graph(); g.add_nodes_from(range(len(self.H)))
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

# ---------------------------------------------------------------------------
def opt_cpsat(inst, respect=None, free=None, forbid_edit=None, time_limit=60, upper=None, member_ze=False):
    """minimum number of edits of a cluster-editing solution.
    respect: set of pairs that must keep their G-status (unedited); free: pairs that cost nothing (and are unconstrained apart from transitivity);
    returns (cost, partition-as-pair-status dict) or None if infeasible."""
    n = inst.n
    m = cp_model.CpModel()
    x = {}
    for u, v in itertools.combinations(range(n), 2):
        x[(u, v)] = m.NewBoolVar(f"x{u}_{v}")   # 1 = together
    for a, b, c in itertools.combinations(range(n), 3):
        m.AddBoolOr([x[(a, b)].Not(), x[(b, c)].Not(), x[(a, c)]])
        m.AddBoolOr([x[(a, b)].Not(), x[(a, c)].Not(), x[(b, c)]])
        m.AddBoolOr([x[(a, c)].Not(), x[(b, c)].Not(), x[(a, b)]])
    cost = []
    respect = set(respect or ())
    free = set(free or ())
    for p in x:
        e = p in inst.E
        if p in respect:
            m.Add(x[p] == (1 if e else 0))
        elif p in free:
            pass
        else:
            cost.append(x[p].Not() if e else x[p])
    if member_ze:
        for (a, b, c) in inst.H:
            m.AddBoolOr([x[key(a, b)], x[key(b, c)]])
    m.Minimize(sum(cost))
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = time_limit
    s.parameters.num_workers = 4
    st = s.Solve(m)
    if st == cp_model.INFEASIBLE:
        return None
    if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        raise RuntimeError("cp-sat status %s" % s.StatusName(st))
    if st != cp_model.OPTIMAL:
        raise RuntimeError("not optimal")
    return int(round(s.ObjectiveValue())), {p: s.Value(x[p]) for p in x}

def ze_states(inst, respect_hard=True):
    """all ZE state vectors (tuples over {0,1,2}) of the member set: every hard pair keeps its status; DFS with triple tests."""
    n = inst.n; H = inst.H; m = len(H)
    # triples that involve at least one soft pair; all-hard triples are static
    soft_t = []
    for x_, y_, z_ in itertools.combinations(range(n), 3):
        ps = [key(x_, y_), key(y_, z_), key(x_, z_)]
        owners = [inst.owner.get(p) for p in ps]
        own = sorted({o[0] for o in owners if o})
        if not own:
            if sum(inst.adj[p[0]][p[1]] for p in ps) == 2:
                return []          # static violation
            continue
        soft_t.append((ps, owners, own))
    order = list(range(m))
    pos = {i: k for k, i in enumerate(order)}
    by_last = [[] for _ in range(m)]
    for ps, owners, own in soft_t:
        by_last[max(pos[i] for i in own)].append((ps, owners))
    X = {0: (0, 1, 1), 1: (1, 0, 1), 2: (0, 0, 1)}   # togetherness of (ab,bc,ac) per state
    res = []
    st = [None] * m
    def ok(k):
        for ps, owners in by_last[k]:
            s = 0
            for p, o in zip(ps, owners):
                if o is None:
                    s += inst.adj[p[0]][p[1]]
                else:
                    s += X[st[o[0]]][o[1]]
            if s == 2:
                return False
        return True
    def rec(k):
        if k == m:
            res.append(tuple(st)); return
        for s in (0, 1, 2):
            st[k] = s
            if ok(k):
                rec(k + 1)
        st[k] = None
    rec(0)
    return res

def edit_cost_of_states(inst, states):
    """cost of a ZE state vector = |H| (sanity)"""
    return len(inst.H)


def opt_csp_cpsat(inst, T=None, time_limit=60):
    """csp-variant relaxation F'(T): min over assignments of statuses to soft+local pairs (far pairs in T fixed at their G status,
    far pairs outside T absent and every triple containing one is dropped) such that all remaining triples are not induced P3.
    T=None: all far pairs respected.  Independent of w9b_leaf (CP-SAT instead of bucket elimination)."""
    n = inst.n
    far = set(inst.far_pairs())
    respected = far if T is None else set(T)
    m = cp_model.CpModel()
    x = {}
    cost = []
    for u, v in itertools.combinations(range(n), 2):
        p = (u, v)
        if p in far:
            continue
        x[p] = m.NewBoolVar(f"x{u}_{v}")
        e = p in inst.E
        # togetherness x; edit cost if x != e ; soft pairs: cost counted too (members' edits)
        cost.append(x[p].Not() if e else x[p])
    def status(p):
        """returns ('var',lit) or ('const',0/1) or None (dropped)"""
        if p in x: return ('var', x[p])
        if p in respected: return ('const', 1 if p in inst.E else 0)
        return None
    for a, b, c in itertools.combinations(range(n), 3):
        ps = [key(a, b), key(b, c), key(a, c)]
        st = [status(p) for p in ps]
        if any(s is None for s in st): continue
        # not exactly two together
        lits = []; const = 0
        for s in st:
            if s[0] == 'const': const += s[1]
            else: lits.append(s[1])
        k = len(lits)
        # sum(lits)+const != 2
        if k == 0:
            if const == 2:
                return None
            continue
        t = 2 - const
        if t < 0 or t > k: continue
        # forbid sum(lits) == t : encode via clauses over assignments (k<=3)
        for assign in itertools.product((0, 1), repeat=k):
            if sum(assign) == t:
                m.AddBoolOr([(l.Not() if v else l) for l, v in zip(lits, assign)])
    m.Minimize(sum(cost))
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = time_limit; s.parameters.num_workers = 4
    st = s.Solve(m)
    if st == cp_model.INFEASIBLE: return None
    if st != cp_model.OPTIMAL: raise RuntimeError("cp-sat %s" % s.StatusName(st))
    return int(round(s.ObjectiveValue()))

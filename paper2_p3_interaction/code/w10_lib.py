"""Wave 10 library: thick one-hot chains OH_m + guest calculus (permanent guests, lock guests, comparator guests, tunnels).

Thin chain OH_m: code/w8c_oh.py (ZE solutions: spine word R^a M L^(m-1-a), a = 0..m-1 = 'position').
Thick: every thin vertex v -> super-vertex (tag, v, i), i<p, glued by a hard clique minus pad pairs; every thin member -> p^2 Latin members.
A *guest* is a single new vertex with a set W of super-vertices it is soft to (pad members: one per pair of a 1-factor of the super-vertex).
Semantics (all VERIFIED-BY-COMPUTER here, argued in notes/wave10_hardness_proof.md):
  guest joins a base cluster C  <=>  every super-vertex of C is in W  (otherwise hard non-edge inside the cluster).
"""
import sys, os, itertools, functools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w8c_lib import *
from w8c_thick import Builder
from w8c_oh import oh_instance
from w8c_chainbuild import one_factorization, cluster_sets
import networkx as nx


@functools.lru_cache(None)
def thin(m):
    n0, E0, mem0, names0, idx0 = oh_instance(m)
    cs = cluster_sets(m)
    # profile[v][a] = frozenset of thin names in the cluster of v at position a, minus v
    prof = {v: {} for v in names0}
    for a, clusters in cs.items():
        for c in clusters:
            for v in c:
                prof[v][a] = frozenset(c - {v})
    return n0, E0, mem0, names0, idx0, prof


class World:
    def __init__(self, p, thin_mode=False):
        self.p = p
        self.thin_mode = thin_mode
        self.free = set()          # thin mode: free (guest, super-vertex) pairs as name pairs
        self.B = Builder()
        self.fact = one_factorization(p)
        self.used = {}
        self.pads = {}
        self.chains = {}
        self.sv_names = set()

    def take_matching(self, sv):
        k = self.used.get(sv, 0)
        assert k < len(self.fact), ('out of matchings at', sv)
        self.used[sv] = k + 1
        M = self.fact[k]
        self.pads.setdefault(sv, set()).update(M)
        return M

    def add_chain(self, tag, m):
        ch = Chain(self, tag, m)
        self.chains[tag] = ch
        return ch

    def guest_pads(self, y, svs):
        if self.thin_mode:
            for (tag, v) in svs:
                self.free.add(frozenset([y, (tag, v)]))
            return
        for sv in sorted(svs, key=str):
            tag, v = sv
            for (i, k) in self.take_matching(sv):
                self.B.member(y, (tag, v, i), (tag, v, k))

    def finish(self):
        p = self.p
        if self.thin_mode:
            return self.B.finish()
        for sv in self.sv_names:
            tag, v = sv
            padset = self.pads.get(sv, set())
            g = nx.Graph(); g.add_nodes_from(range(p))
            for i in range(p):
                for k in range(i + 1, p):
                    if (i, k) not in padset:
                        g.add_edge(i, k); self.B.hard((tag, v, i), (tag, v, k))
            assert nx.is_connected(g), ('glue disconnected', sv, len(padset))
        return self.B.finish()


class Chain:
    def __init__(self, W, tag, m):
        self.W, self.tag, self.m = W, tag, m
        self.n0, self.E0, self.mem0, self.names0, self.idx0, self.prof = thin(m)
        p = W.p
        soft0 = set()
        for t in self.mem0:
            soft0 |= pairs_of(t)
        if W.thin_mode:
            for (a, b, c) in self.mem0:
                W.B.member((tag, self.names0[a]), (tag, self.names0[b]), (tag, self.names0[c]))
            for (u, v) in self.E0:
                if (u, v) not in soft0:
                    W.B.hard((tag, self.names0[u]), (tag, self.names0[v]))
            return
        for (a, b, c) in self.mem0:
            A, Bn, C = self.names0[a], self.names0[b], self.names0[c]
            for i in range(p):
                for j in range(p):
                    W.B.member((tag, A, i), (tag, Bn, j), (tag, C, (i + j) % p))
        for (u, v) in self.E0:
            if (u, v) in soft0:
                continue
            for i in range(p):
                for j in range(p):
                    W.B.hard((tag, self.names0[u], i), (tag, self.names0[v], j))
        for v in self.names0:
            W.sv_names.add((tag, v))

    # ---- profile calculus
    def allowed(self, host, forbid):
        """thin names W (as set) such that a guest soft exactly to super-vertices W (plus host) can join the host cluster at
        position a iff a not in the returned effective forbidden set."""
        Wset = {host}
        for a, c in self.prof[host].items():
            if a not in forbid:
                Wset |= set(c)
        eff = {a for a, c in self.prof[host].items() if not set(c) <= Wset}
        return Wset, eff

    def sv(self, thins):
        return {(self.tag, v) for v in thins}

    def perm_guest(self, z, host, forbid):
        """z forced into the host cluster; soft to everything the host cluster may contain at positions not in `forbid`.
        Returns the effective forbidden set (positions whose host cluster contains a vertex z is not soft to)."""
        W, tag, B = self.W, self.tag, self.W.B
        if W.thin_mode:
            B.hard((tag, host), z)
        else:
            for (i, k) in W.take_matching((tag, host)):
                B.member((tag, host, i), z, (tag, host, k))
        Wset, eff = self.allowed(host, forbid)
        W.guest_pads(z, self.sv(Wset - {host}))
        return eff

    def pos_vertices(self, a):
        if self.W.thin_mode:
            return [(self.tag, ('s', 2 * a)), (self.tag, ('s', 2 * a + 1)), (self.tag, ('s', 2 * a + 2))]
        return [(self.tag, ('s', 2 * a), 0), (self.tag, ('s', 2 * a + 1), 0), (self.tag, ('s', 2 * a + 2), 0)]


def lock_guest(ch, y, x, z, host, forbid):
    """y closable guest hosted at `host`; member (x,y,z); anchor x.  Locking x (x singleton / tunnel respected) forces y~z, i.e. y in host cluster."""
    Wset, eff = ch.allowed(host, forbid)
    ch.W.guest_pads(y, ch.sv(Wset))
    ch.W.B.member(x, y, z)
    return eff


def cmp_guest(chA, hostA, forbA, chB, hostB, forbB, y, zA, zB):
    """comparator guest y: member (zA,y,zB); y joins the A host cluster at positions not in effA, the B host cluster at positions not in effB."""
    WA, effA = chA.allowed(hostA, forbA)
    WB, effB = chB.allowed(hostB, forbB)
    svs = chA.sv(WA) | chB.sv(WB)
    chA.W.guest_pads(y, svs)
    chA.W.B.member(zA, y, zB)
    return effA, effB


class Inst:
    def __init__(self, W):
        self.W = W
        self.n, self.E, self.members = W.finish()
        self.idx = W.B.idx
        self.valid = check_instance(self.n, self.E, self.members)
        self.c = max(sum(1 for t in self.members if v in t) for v in range(self.n))

    def pos_of(self, lab, ch):
        idx = self.idx
        res = []
        for q in range(ch.m):
            a, b, c = [idx[v] for v in ch.pos_vertices(q)]
            if lab[a] == lab[b] == lab[c]:
                res.append(q)
        return res

    def solve(self, singletons=(), together=(), apart=(), drop_hard=(), **kw):
        """ZE solve. singletons/together/apart given as vertex names. drop_hard: pairs of names whose hard edge is removed (tunnel cut)."""
        idx = self.idx
        E = self.E
        if drop_hard:
            E = set(E)
            for (u, v) in drop_hard:
                E.discard(key(idx[u], idx[v]))
        return ze_solve_fast(self.n, E, self.members,
                             singletons=[idx[x] for x in singletons],
                             fixed_together=[(idx[u], idx[v]) for u, v in together],
                             fixed_apart=[(idx[u], idx[v]) for u, v in apart], **kw)

    def force_pos(self, ch, a):
        v = ch.pos_vertices(a)
        return [(v[0], v[1]), (v[1], v[2])]


def ze_atm(n, E, members, free, singletons=(), fixed_apart=(), fixed_together=(), all_solutions=False, limit=None):
    """ZE decision in the abstract thin-guest model: vertices 0..n-1, members (one edit each), `free` pairs (vertex index pairs, no cost, no constraint
    other than transitivity), all other non-member pairs hard (edge iff in E).  Returns labels or None."""
    soft = set()
    for m in members:
        soft |= pairs_of(m)
    vp = VarPool()
    var = {}
    for pq in soft | free:
        var[pq] = vp.new()
    TRUE = vp.new()
    cl = [[TRUE]]
    def lit(u, v):
        pq = key(u, v)
        if pq in var:
            return var[pq]
        return TRUE if pq in E else -TRUE
    adj = [set() for _ in range(n)]
    for pq in var:
        adj[pq[0]].add(pq[1]); adj[pq[1]].add(pq[0])
    for pq in E:
        if pq not in var:
            adj[pq[0]].add(pq[1]); adj[pq[1]].add(pq[0])
    seen = set()
    for v in range(n):
        nb = sorted(adj[v])
        for i in range(len(nb)):
            for j in range(i + 1, len(nb)):
                t = tuple(sorted((nb[i], v, nb[j])))
                if t in seen:
                    continue
                seen.add(t)
                u, v_, w = t
                a, b, c = lit(u, v_), lit(v_, w), lit(u, w)
                if all(abs(x) == TRUE for x in (a, b, c)):
                    if (a == TRUE) + (b == TRUE) + (c == TRUE) == 2:
                        cl.append([-TRUE])
                    continue
                cl += [[-a, -b, c], [-a, -c, b], [-b, -c, a]]
    for (a, b, c) in members:
        e1, e2, e3 = -var[key(a, b)], -var[key(b, c)], var[key(a, c)]
        cl += [[e1, e2, e3], [-e1, -e2], [-e1, -e3], [-e2, -e3]]
    def fix(u, v, val):
        pq = key(u, v)
        if pq in var:
            cl.append([var[pq]] if val else [-var[pq]])
        elif (pq in E) != val:
            cl.append([-TRUE])
    for x in singletons:
        for v in range(n):
            if v != x:
                fix(x, v, False)
    for (u, v) in fixed_apart:
        fix(u, v, False)
    for (u, v) in fixed_together:
        fix(u, v, True)
    s = Cadical153(bootstrap_with=cl)
    out = []
    def decode(model):
        ms = set(l for l in model if l > 0)
        parent = list(range(n))
        def f(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]; x = parent[x]
            return x
        for pq in E:
            if pq not in var:
                parent[f(pq[0])] = f(pq[1])
        for pq, vv in var.items():
            if vv in ms:
                parent[f(pq[0])] = f(pq[1])
        labs = {}
        lab = []
        for v in range(n):
            r = f(v)
            if r not in labs:
                labs[r] = len(labs)
            lab.append(labs[r])
        return lab
    while s.solve():
        model = s.get_model()
        if not all_solutions:
            s.delete()
            return decode(model)
        out.append(decode(model))
        if limit and len(out) >= limit:
            break
        ms = set(l for l in model if l > 0)
        s.add_clause([(-v if v in ms else v) for v in var.values()])
    s.delete()
    return out if all_solutions else None


class ThinInst:
    """abstract thin-guest instance (World(thin_mode=True))"""
    def __init__(self, W):
        self.W = W
        self.n, self.E, self.members = W.finish()
        self.idx = W.B.idx
        self.free = {key(self.idx[a], self.idx[b]) for a, b in (tuple(f) for f in W.free)}
        self.valid = check_instance(self.n, self.E, self.members)

    def solve(self, singletons=(), together=(), apart=(), drop_hard=(), **kw):
        idx = self.idx
        E = self.E
        if drop_hard:
            E = set(E)
            for (u, v) in drop_hard:
                E.discard(key(idx[u], idx[v]))
        return ze_atm(self.n, E, self.members, self.free,
                      singletons=[idx[x] for x in singletons],
                      fixed_together=[(idx[u], idx[v]) for u, v in together],
                      fixed_apart=[(idx[u], idx[v]) for u, v in apart], **kw)

    def force_pos(self, ch, a):
        v = ch.pos_vertices(a)
        return [(v[0], v[1]), (v[1], v[2])]

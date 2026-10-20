"""w9b_leaf.py -- the leaf oracle 'no far edit' as a CSP over bag-local variables, solved by bucket elimination
(variable elimination along a min-fill order, which realises a tree decomposition; the induced width is reported).

Variables:  ('m',i)   member i, domain = the 7 non-empty subsets of {ab,bc,ac} (bitmask 1..7), cost popcount
            ('b',u,v) local hard pair, domain {0,1} (1 = edited), cost 1
Constraints: every vertex triple is not an induced P3 of the edited graph; a far hard pair is NEVER edited (static),
             far pairs outside T are dropped together with every triple containing them (csp variant, option T)
leaf_cost(inst, T=None): min total cost, or None.  T=None means all far pairs respected (the leaf test F(N)).
"""
import itertools
import numpy as np
from w9b_lib import Inst, key

INF = 10 ** 6
SUBS = list(range(1, 8))            # bit0 = ab edited, bit1 = bc edited, bit2 = ac edited


class CSP:
    def __init__(self, inst, T=None):
        self.inst = inst
        far = set(inst.far_pairs())
        self.far = far
        respected = far if T is None else set(T)
        self.vars = []                # list of ('m',i) / ('b',u,v)
        self.dom = {}
        self.cost = {}
        for i in range(len(inst.H)):
            v = ('m', i); self.vars.append(v); self.dom[v] = 7; self.cost[v] = np.array([bin(s).count('1') for s in SUBS], float)
        self.locbits = {}
        for p in inst.local_hard_pairs():
            v = ('b',) + p; self.vars.append(v); self.dom[v] = 2; self.cost[v] = np.array([0., 1.]); self.locbits[p] = v
        self.factors = {}             # scope(tuple of vars, sorted) -> ndarray (additive, 0 / INF)
        n = inst.n
        self.dropped = 0
        for x, y, z in itertools.combinations(range(n), 3):
            ps = [key(x, y), key(y, z), key(x, z)]
            terms = []          # each: ('s', var, bit, base) | ('l', var, base) | ('c', base)
            skip = False
            for p in ps:
                e = 1 if inst.adj[p[0]][p[1]] else 0
                o = inst.owner.get(p)
                if o is not None:
                    terms.append(('s', ('m', o[0]), o[1], e))
                elif p in far:
                    if p in respected:
                        terms.append(('c', e))
                    else:
                        skip = True; break
                else:
                    terms.append(('l', self.locbits[p], e))
            if skip:
                self.dropped += 1
                continue
            scope = sorted({t[1] for t in terms if t[0] != 'c'})
            if not scope:
                if sum(t[1] for t in terms) == 2:
                    self.factors[()] = np.array(INF)       # static violation
                continue
            shape = tuple(self.dom[v] for v in scope)
            tot = np.zeros(shape, int)
            for t in terms:
                if t[0] == 'c':
                    tot += t[1]
                elif t[0] == 's':
                    _, v, bit, e = t
                    ax = scope.index(v)
                    arr = np.array([e ^ ((s >> bit) & 1) for s in SUBS])
                    sh = [1] * len(scope); sh[ax] = 7
                    tot = tot + arr.reshape(sh)
                else:
                    _, v, e = t
                    ax = scope.index(v)
                    arr = np.array([e ^ 0, e ^ 1])
                    sh = [1] * len(scope); sh[ax] = 2
                    tot = tot + arr.reshape(sh)
            f = np.where(tot == 2, INF, 0).astype(float)
            sc = tuple(scope)
            if sc in self.factors: self.factors[sc] = self.factors[sc] + f
            else: self.factors[sc] = f

    def solve(self):
        # factors list incl. unary costs
        facs = [(sc, f) for sc, f in self.factors.items() if sc != ()]
        if () in self.factors and self.factors[()] >= INF:
            return None, 0
        for v in self.vars:
            facs.append(((v,), self.cost[v].copy()))
        # min-fill elimination order
        adj = {v: set() for v in self.vars}
        for sc, f in facs:
            for a in sc:
                for b in sc:
                    if a != b: adj[a].add(b)
        remaining = set(self.vars)
        width = 0
        const = 0.0
        while remaining:
            best = None
            for v in remaining:
                nb = adj[v] & remaining
                fill = sum(1 for a, b in itertools.combinations(nb, 2) if b not in adj[a])
                sz = 1
                for u in nb: sz *= self.dom[u]
                key_ = (fill, len(nb))
                if best is None or key_ < best[0]: best = (key_, v)
            v = best[1]
            nb = sorted(adj[v] & remaining)
            width = max(width, len(nb))
            rel = [(sc, f) for sc, f in facs if v in sc]
            facs = [(sc, f) for sc, f in facs if v not in sc]
            # combine
            scope = sorted({u for sc, f in rel for u in sc})
            shape = tuple(self.dom[u] for u in scope)
            tot = np.zeros(shape)
            for sc, f in rel:
                sh = [1] * len(scope)
                for u in sc: sh[scope.index(u)] = self.dom[u]
                tot = tot + f.reshape(sh)
            ax = scope.index(v)
            red = tot.min(axis=ax)
            newscope = tuple(u for u in scope if u != v)
            if newscope:
                facs.append((newscope, red))
            else:
                const += float(red)
            remaining.discard(v)
            for a in nb:
                adj[a].discard(v)
                for b in nb:
                    if a != b: adj[a].add(b)
        if const >= INF: return None, width
        return int(round(const)), width


def leaf_cost(inst, T=None):
    c = CSP(inst, T)
    return c.solve()


def min_certificate(inst, budget, N=None, order=None):
    """smallest-by-deletion minimal certificate T subset N with csp-variant F'(T) > budget (assumes F'(N) > budget)."""
    far = list(inst.far_pairs()) if N is None else list(N)
    T = list(far)
    for q in list(T):
        T2 = [p for p in T if p != q]
        c, _ = leaf_cost(inst, T2)
        if c is None or c > budget:
            T = T2
    return T

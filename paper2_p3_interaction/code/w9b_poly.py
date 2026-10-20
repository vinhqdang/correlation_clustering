"""w9b_poly.py -- polymorphisms of the zero-excess constraint language (wave 9b, written from scratch).

Domain D = {0,1,2} = {L,R,M} (member (a,b,c): L = edit ab (b~c, a apart), R = edit bc (a~b), M = edit ac (all together)).
Indicator of 'pair together' after the single edit:  X[e][s],  e in {ab,bc,ac}.
Every triple of vertices gives 'not exactly two of its three pairs together' (transitivity).  A term is a pair of a member
(variable + role e) or a constant 0/1 (hard pair).  Relations:
  unary lists, binary NAND/EQ relations (eps=0/1), ternary soft-triangle relations, and conjunctions of single-triple relations.
Polymorphism search: SAT with CEGAR (counterexample guided clause generation).  pysat only.
"""
import itertools, sys, time
from pysat.solvers import Cadical153

D = (0, 1, 2)
X = {'ab': (0, 1, 1), 'bc': (1, 0, 1), 'ac': (0, 0, 1)}   # X[e][state]
ROLES = ('ab', 'bc', 'ac')

def triple_rel(terms, nv):
    """terms: list of ('v',i,role) or ('c',val); nv variables. relation = set of tuples in D^nv with sum != 2"""
    out = set()
    for t in itertools.product(D, repeat=nv):
        s = 0
        for tm in terms:
            s += X[tm[2]][t[tm[1]]] if tm[0] == 'v' else tm[1]
        if s != 2:
            out.add(t)
    return frozenset(out)

def single_triples():
    """all single-triple relations, as dict name -> (arity, relation)"""
    rels = {}
    for e in ROLES:                       # unary: two constants
        for c in (0, 1, 2):
            rels[('u', e, c)] = (1, triple_rel([('v', 0, e), ('c', c)] if False else [('v', 0, e)] + [('c', 1)] * c, 1))
    for e1 in ROLES:
        for e2 in ROLES:
            for eps in (0, 1):
                rels[('b', e1, e2, eps)] = (2, triple_rel([('v', 0, e1), ('v', 1, e2), ('c', eps)], 2))
    for e1, e2, e3 in itertools.product(ROLES, repeat=3):
        rels[('t', e1, e2, e3)] = (3, triple_rel([('v', 0, e1), ('v', 1, e2), ('v', 2, e3)], 3))
    return rels

def lists():
    rels = {}
    for k in (1, 2):
        for S in itertools.combinations(D, k):
            rels[('L', S)] = (1, frozenset((s,) for s in S))
    return rels

# ----------------------------------------------------------------------------
class PolySearch:
    """find f: D^r -> D preserving all relations, with extra constraints."""
    def __init__(self, r, rels, nu=False, idem=True, conservative=False, extra=None):
        self.r = r
        self.tuples = list(itertools.product(D, repeat=r))
        self.tid = {t: i for i, t in enumerate(self.tuples)}
        self.nv = len(self.tuples) * 3
        self.rels = rels
        self.S = Cadical153()
        # one-hot
        for i in range(len(self.tuples)):
            self.S.add_clause([self.var(i, v) for v in D])
            for v, w in itertools.combinations(D, 2):
                self.S.add_clause([-self.var(i, v), -self.var(i, w)])
        if idem:
            for v in D:
                self.S.add_clause([self.var(self.tid[(v,) * r], v)])
        if conservative:
            for t in self.tuples:
                allowed = set(t)
                for v in D:
                    if v not in allowed:
                        self.S.add_clause([-self.var(self.tid[t], v)])
        if nu:
            for t in self.tuples:
                from collections import Counter
                c = Counter(t)
                v, m = c.most_common(1)[0]
                if m >= r - 1:
                    self.S.add_clause([self.var(self.tid[t], v)])
        if extra:
            extra(self)
        self.added = 0

    def var(self, ti, v):
        return ti * 3 + v + 1

    def val(self, model, ti):
        for v in D:
            if model[self.var(ti, v) - 1] > 0:
                return v

    def check(self, f):
        """f: list of values by tuple index. return violating (rel name, rows) or None. Exhaustive over rows."""
        r = self.r
        for name, (k, R) in self.rels.items():
            Rl = list(R)
            for rows in itertools.product(Rl, repeat=r):
                img = tuple(f[self.tid[tuple(rows[j][i] for j in range(r))]] for i in range(k))
                if img not in R:
                    return name, rows
        return None

    def add_violation(self, name, rows):
        k, R = self.rels[name]
        r = self.r
        cols = [self.tid[tuple(rows[j][i] for j in range(r))] for i in range(k)]
        for bad in itertools.product(D, repeat=k):
            if bad not in R:
                self.S.add_clause([-self.var(cols[i], bad[i]) for i in range(k)])
        self.added += 1

    def solve(self, maxit=100000, verbose=False):
        it = 0
        while it < maxit:
            it += 1
            if not self.S.solve():
                return None
            m = self.S.get_model()
            f = [self.val(m, i) for i in range(len(self.tuples))]
            v = self.check(f)
            if v is None:
                return f
            # add several violations at once: all violating rows for this relation is expensive; add the found one
            self.add_violation(*v)
            if verbose and it % 200 == 0:
                print('it', it, 'clauses added', self.added, flush=True)
        return 'timeout'

def describe(f, r):
    return {t: v for t, v in zip(itertools.product(D, repeat=r), f)}

if __name__ == '__main__':
    r = int(sys.argv[1]); which = sys.argv[2]
    rels = {}
    rels.update(lists())
    st = single_triples()
    if which == 'binary':
        rels.update({k: v for k, v in st.items() if v[0] <= 2})
    else:
        rels.update(st)
    t0 = time.time()
    ps = PolySearch(r, rels, nu=True, conservative=True)
    res = ps.solve(verbose=True)
    print('r', r, which, 'result', 'NONE' if res is None else ('timeout' if res == 'timeout' else 'FOUND'), 'clauses', ps.added, 'time', round(time.time() - t0, 1))

# ---------------------------------------------------------------------------
# realisable language: scope relations at a shared vertex (partial conjunctions) and soft-triangle ternary relations
ROLEPAIRS = [('ab', 'ac'), ('ab', 'bc'), ('bc', 'ac')]   # roles of the two pairs {q,p1},{q,p2} for q at position a / b / c  (q=a:ab,ac ; q=b:ab,bc ; q=c: bc,ac)

def cross_relation(rA, rB, eps, S):
    """binary relation on (sA,sB): conjunction over (i,j) in S of  x_rA[i](sA)+x_rB[j](sB)+eps[i][j] != 2"""
    out = set()
    for sA in D:
        for sB in D:
            ok = True
            for (i, j) in S:
                if X[rA[i]][sA] + X[rB[j]][sB] + eps[i][j] == 2:
                    ok = False; break
            if ok: out.add((sA, sB))
    return frozenset(out)

def binary_language(max_soft_pattern='all'):
    """returns dict name -> (2, rel). pattern: 'none' (full conj only), 'c2' (soft cross pairs form a matching), 'all'"""
    cross = [(0, 0), (0, 1), (1, 0), (1, 1)]
    matchings = [(), ((0, 0),), ((0, 1),), ((1, 0),), ((1, 1),), ((0, 0), (1, 1)), ((0, 1), (1, 0))]
    out = {}
    allsub = [tuple(c) for k in range(5) for c in itertools.combinations(cross, k)]
    for rA in ROLEPAIRS:
        for rB in ROLEPAIRS:
            for eps in itertools.product((0, 1), repeat=4):
                e = [[eps[0], eps[1]], [eps[2], eps[3]]]
                if max_soft_pattern == 'none':
                    softs = [()]
                elif max_soft_pattern == 'c2':
                    softs = matchings
                else:
                    softs = allsub
                for sp in softs:
                    S = tuple(x for x in cross if x not in sp)
                    R = cross_relation(rA, rB, e, S)
                    out[('B', rA, rB, eps, sp)] = (2, R)
    # dedupe by relation
    ded = {}
    for k, (a, R) in out.items():
        ded.setdefault(R, k)
    return {k: (2, R) for R, k in ded.items()}

def ternary_language():
    out = {}
    for e1, e2, e3 in itertools.product(ROLES, repeat=3):
        out[('T', e1, e2, e3)] = (3, triple_rel([('v', 0, e1), ('v', 1, e2), ('v', 2, e3)], 3))
    return out

def language(name):
    rels = {}
    rels.update(lists())
    if name == 'free':
        rels.update(single_triples())
    elif name == 'sf':
        rels.update(binary_language('none'))
    elif name == 'c2':
        rels.update(binary_language('c2')); rels.update(ternary_language())
    elif name == 'c3':
        rels.update(binary_language('all')); rels.update(ternary_language())
    elif name == 'ternary':
        rels.update(ternary_language())
    return rels

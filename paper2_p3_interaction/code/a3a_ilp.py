"""ILP (scipy/HiGHS) audit of the same statements for n = 10..13 on adversarial families.
usage: python3 a3a_ilp.py seed ninst nmin nmax [hyp]"""
import sys, random, itertools, collections, time
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix
import a3a_lib as L

_tri = {}
def trans(n):
    if n in _tri: return _tri[n]
    pi = L.pidx(n); rows = []
    for i, j, k in itertools.combinations(range(n), 3):
        a, b, c = pi[(i, j)], pi[(j, k)], pi[(i, k)]
        rows += [(a, b, c), (a, c, b), (b, c, a)]   # z_a + z_b - z_c <= 1
    M = lil_matrix((len(rows), n * (n - 1) // 2))
    for r, (a, b, c) in enumerate(rows):
        M[r, a] = 1; M[r, b] = 1; M[r, c] = -1
    _tri[n] = M.tocsr(); return _tri[n]

class Model:
    def __init__(self, A):
        self.A = A; n = self.n = len(A); self.pi = L.pidx(n); m = n * (n - 1) // 2
        self.ed = np.array([A[i, j] for i, j in itertools.combinations(range(n), 2)], float)
        self.cvec = 1 - 2 * self.ed        # cost = sum(ed) + cvec . z
        self.const = self.ed.sum(); self.T = trans(n); self.m = m
        self.opt = None
    def solve(self, obj, fixed=None, cost_le=None, extra=None):
        """minimise obj.z ; fixed = {pair index: 0/1}; cost_le = max total cost; extra = list of (vec, lo, hi)"""
        lb = np.zeros(self.m); ub = np.ones(self.m)
        for k, v in (fixed or {}).items(): lb[k] = ub[k] = v
        cons = [LinearConstraint(self.T, -np.inf, 1)]
        if cost_le is not None: cons.append(LinearConstraint(self.cvec.reshape(1, -1), -np.inf, cost_le - self.const))
        for vec, lo, hi in (extra or []): cons.append(LinearConstraint(np.array(vec, float).reshape(1, -1), lo, hi))
        r = milp(obj, constraints=cons, integrality=np.ones(self.m), bounds=Bounds(lb, ub))
        if r.status != 0: return None
        return r
    def OPT(self):
        if self.opt is None:
            r = self.solve(self.cvec); self.opt = int(round(r.fun + self.const))
        return self.opt

def gen_clusters(rng, n):
    """2-3 noisy clusters + attached twin classes (adversarial: many classes, shared hubs, equal attachments)"""
    A = L.sym(n)
    k = rng.randint(max(4, n // 2 - 1), n - 2)
    cut = sorted(rng.sample(range(1, k), rng.randint(1, 2)))
    grp = []; prev = 0
    for c in cut + [k]: grp.append(list(range(prev, c))); prev = c
    for g in grp:
        for i, j in itertools.combinations(g, 2):
            if rng.random() < .85: L.add(A, i, j)
    for _ in range(rng.randint(0, 3)):
        i, j = rng.sample(range(k), 2); L.add(A, i, j, 1 - A[i, j])
    v = k; classes = []
    while v < n:
        s = min(rng.randint(1, 4), n - v); cl = list(range(v, v + s)); v += s; classes.append(cl)
        for a, b in itertools.combinations(cl, 2): L.add(A, a, b)
    mode = rng.choice(['disjoint', 'hub', 'same', 'any'])
    hub = rng.randrange(k); same = rng.sample(range(k), rng.randint(1, 2)); used = set()
    for cl in classes:
        if mode == 'disjoint':
            av = [z for z in range(k) if z not in used]; att = rng.sample(av, min(len(av), rng.randint(1, 2))) if av else []
            used |= set(att)
        elif mode == 'hub': att = [hub]
        elif mode == 'same': att = same
        else: att = rng.sample(range(k), rng.randint(1, 3))
        for z in att:
            for x in cl: L.add(A, x, z)
    return A, k

def gen_packed(rng, n):
    """plant a pair-disjoint packing, attach twin classes to its vertices, add edges among V(H); the greedy
    completion (start = planted P3s) keeps the planted R wherever the pair-maximality allows."""
    m = rng.randint(2, 4); share = rng.random() < .4
    verts = 0; P = []
    for i in range(m):
        if share and P and rng.random() < .6:
            a = rng.choice([v for t in P for v in t]); b = verts; c = verts + 1; verts += 2
            t = (b, a, c) if rng.random() < .5 else (a, b, c)   # a is an end or the centre
            P.append(t if True else t)
        else:
            P.append((verts, verts + 1, verts + 2)); verts += 3
    if verts >= n - 1: return None
    A = L.sym(n)
    for (x, y, z) in P:
        L.add(A, x, y); L.add(A, y, z)
    VHs = sorted(set(v for t in P for v in t))
    v = verts; classes = []
    while v < n:
        s = min(rng.randint(1, 4), n - v); cl = list(range(v, v + s)); v += s; classes.append(cl)
        for a, b in itertools.combinations(cl, 2): L.add(A, a, b)
    mode = rng.choice(['disjoint', 'disjoint', 'hub', 'any']); hub = rng.choice(VHs); used = set()
    for cl in classes:
        if mode == 'disjoint':
            av = [z for z in VHs if z not in used]; att = rng.sample(av, min(len(av), rng.randint(1, 2))) if av else []
            used |= set(att)
        elif mode == 'hub': att = [hub]
        else: att = rng.sample(VHs, rng.randint(1, 2))
        for z in att:
            for x in cl: L.add(A, x, z)
    p = rng.choice([0, .1, .25])
    for i, j in itertools.combinations(VHs, 2):
        if rng.random() < p and not any({i, j} <= set(t) for t in P): L.add(A, i, j, 1 - A[i, j])
    keep = []
    for t in P:   # keep planted P3s that are still induced and pair-disjoint
        if A[t[0], t[1]] and A[t[1], t[2]] and not A[t[0], t[2]] if False else None: pass
    return A, P

def instance(rng, n):
    r = rng.random()
    if r < .5:
        g = gen_packed(rng, n)
        if g is not None:
            A, P = g
            used = set(); start = []
            for t in P:
                # recompute orientation: find the centre in the current graph
                for c in t:
                    e = [x for x in t if x != c]
                    if A[c, e[0]] and A[c, e[1]] and not A[e[0], e[1]]:
                        tt = (e[0], c, e[1]); ps = L.tri_pairs(tt)
                        if all(p not in used for p in ps): start.append(tt); used.update(ps)
                        break
            H = L.greedy_packing(A, rng, None, start=start)
            return 'packed', A, H
    if r < .75:
        A, k = gen_clusters(rng, n)
        key = (lambda t: sum(1 for x in t if x >= k)) if rng.random() < .6 else None
        return 'clusters', A, L.greedy_packing(A, rng, key)
    f, A, H = L.make_instance(rng, n)
    return f, A, H

def run(seed, ninst, nmin, nmax, hyp='pairmax', tlimit=600):
    rng = random.Random(seed); C = collections.Counter(); fails = []; t0 = time.time()
    for it in range(ninst):
        if time.time() - t0 > tlimit: break
        n = rng.randint(nmin, nmax)
        f, A, H = instance(rng, n)
        if hyp == 'vertexmax': H = L.vertex_greedy_packing(A, rng)
        if hyp == 'nonmax': H = H[:rng.randint(0, max(0, len(H) - 1))]
        VH, R, Q, classes = L.info(A, H)
        if len(R) < 2 or not H: continue
        pm = L.is_pair_maximal(A, H)
        if hyp == 'pairmax': assert pm
        l4 = L.lemma4_holds(A, R)
        M = Model(A); OPT = M.OPT(); lopt = OPT - len(H); pi = M.pi
        C['inst'] += 1; C['pairmax' if pm else 'notpairmax'] += 1
        C['l=%d' % min(lopt, 5)] += 1
        C['R>=4'] += len(R) >= 4; C['classes>=3'] += len(classes) >= 3
        C['multi_vertex_class'] += any(len(X) >= 2 for X in classes)
        if not l4: C['lemma4_fail'] += 1
        # twin lemma on every optimum, all true twin pairs of G
        tw = L.true_twin_pairs(A)
        if tw:
            obj = np.zeros(M.m)
            for u, v in tw: obj[pi[(u, v)]] = 1       # minimise sum z -> maximise separations
            r = M.solve(obj, cost_le=OPT)
            sep = len(tw) - int(round(r.fun))
            C['twin_pairs'] += len(tw)
            if sep > 0: C['TWIN_FAIL'] += 1; fails.append(('twin', A.tolist(), H))
        # C1 'every optimum' form: max R-R edits among optima ; and 'some optimum' form: min R-R edits
        rr = list(itertools.combinations(R, 2)); obj = np.zeros(M.m); const = 0
        for a, b in rr:
            k = pi[(a, b)]
            if A[a, b]: obj[k] = 1; const += 0      # edge deleted iff z=0 -> edits = 1 - z ; minimise edits => obj -z
            if A[a, b]: obj[k] = -1; const += 1
            else: obj[k] = 1
        rmin = M.solve(obj, cost_le=OPT); some_edits = int(round(rmin.fun + const))
        rmax = M.solve(-obj, cost_le=OPT); every_edits = int(round(-rmax.fun + const))
        if some_edits > 0: C['C1_some_FAIL'] += 1; fails.append(('c1some', hyp, A.tolist(), H))
        if every_edits > 0: C['C1_every_optimum_FAIL'] += 1
        # N3 / P2 / P4 per class
        qmax = max(len(Q[v]) for v in range(n))
        for X in classes:
            AX = [z for z in range(n) if z not in X and A[X[0], z]]
            x0 = X[0]
            if len(X) > lopt:
                obj = np.zeros(M.m)
                for v in range(n):
                    if v == x0 or v in X: continue
                    obj[pi[(x0, v)]] = -1 if A[x0, v] else 1    # edits at x0
                r = M.solve(-obj, cost_le=OPT); e = int(round(-r.fun + sum(1 for v in AX)))
                C['N3_tests'] += 1
                if e > 0: C['N3_FAIL'] += 1; fails.append(('N3', A.tolist(), H, X))
            # P2: among optima with z_{x0,AX}=0 maximise |B|
            outside = [a for a in VH if a not in AX]
            obj = np.zeros(M.m)
            for a in outside: obj[pi[(x0, a)]] = -1
            fixed = {pi[(x0, z)]: 0 for z in AX}
            r = M.solve(obj, fixed=fixed, cost_le=OPT)
            C['P2_tests'] += 1
            if r is not None and -r.fun > 0.5: C['P2_FAIL'] += 1; fails.append(('P2', A.tolist(), H, X))
            if len(AX) >= lopt + qmax:
                r = M.solve(obj, cost_le=OPT); C['P4_tests'] += 1
                if -r.fun > 0.5: C['P4_FAIL'] += 1; fails.append(('P4', A.tolist(), H, X))
            # N5
            if len(X) >= 2:
                x = X[-1]; keep = [v for v in range(n) if v != x]
                O2 = Model(A[np.ix_(keep, keep)]).OPT(); C['N5_tests'] += 1
                if O2 - len(H) <= len(X) - 2:
                    C['N5_applicable'] += 1
                    if O2 != OPT: C['N5_FAIL'] += 1; fails.append(('N5', A.tolist(), H, X))
                elif O2 != OPT: C['OPT_changes(nonapplicable)'] += 1
    return C, fails

if __name__ == '__main__':
    seed, ninst, nmin, nmax = map(int, sys.argv[1:5]); hyp = sys.argv[5] if len(sys.argv) > 5 else 'pairmax'
    C, fails = run(seed, ninst, nmin, nmax, hyp)
    print('seed', seed, 'hyp', hyp, dict(sorted(C.items()))); print('n_fail_records', len(fails))
    for f in fails[:2]: print(f)

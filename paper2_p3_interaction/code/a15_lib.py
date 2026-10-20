"""a15: independent re-check tools for the wave-10 reduction.
 - contract(): contract every super-vertex K_v (p vertices) to one node (valid when no K_v is split, which Lemma 5 guarantees for ex<=l);
   weights e_AB / ne_AB = number of edge / non-edge vertex pairs between classes.
 - wopt(): exact weighted correlation clustering on the contracted instance (CP-SAT, lazy transitivity triples), ALL class pairs editable.
"""
import sys, os, itertools, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w10_red import *
import numpy as np

def classes_of(I):
    names = I.W.B.names
    cls = {}
    cl_id = {}
    for v, nm in enumerate(names):
        if isinstance(nm, tuple) and len(nm) == 3 and isinstance(nm[1], tuple) and nm[1][0] in ('s', 'e', 'cap') or (isinstance(nm, tuple) and len(nm) == 3 and isinstance(nm[2], int) and isinstance(nm[1], tuple) and nm[0] not in ('x', 'y', 'z', 'c')):
            key_ = (nm[0], nm[1])
        else:
            key_ = nm
        if key_ not in cl_id:
            cl_id[key_] = len(cl_id)
        cls[v] = cl_id[key_]
    return cls, cl_id

def contract(I):
    cls, cl_id = classes_of(I)
    N = len(cl_id)
    size = [0] * N
    for v in range(I.n):
        size[cls[v]] += 1
    e = {}
    for (u, v) in I.E:
        a, b = cls[u], cls[v]
        k = (min(a, b), max(a, b))
        e[k] = e.get(k, 0) + 1
    ne = {}
    # non-edges: total pairs - edges
    for a in range(N):
        for b in range(a, N):
            tot = size[a] * (size[a] - 1) // 2 if a == b else size[a] * size[b]
            ne[(a, b)] = tot - e.get((a, b), 0)
    return N, size, e, ne, cl_id

def wopt(N, size, e, ne, bound=None, verbose=True, tlimit=3600, workers=4, hint=None):
    """minimise sum_{A<B} (T? ne : e) + sum_A ne_AA.  returns (cost, labels)"""
    from ortools.sat.python import cp_model
    const = sum(ne[(a, a)] for a in range(N))
    linked = [set() for _ in range(N)]
    for (a, b), c in e.items():
        if a != b and c > 0:
            linked[a].add(b); linked[b].add(a)
    pairs = [(a, b) for a in range(N) for b in range(a + 1, N)]
    M = cp_model.CpModel()
    T = {pq: M.NewBoolVar('t%d_%d' % pq) for pq in pairs}
    def t(a, b):
        return T[(a, b) if a < b else (b, a)]
    seen = set()
    def add_triple(a, b, c):
        s = tuple(sorted((a, b, c)))
        if s in seen: return
        seen.add(s)
        x, y, z = t(s[0], s[1]), t(s[1], s[2]), t(s[0], s[2])
        M.AddBoolOr([x.Not(), y.Not(), z]); M.AddBoolOr([x.Not(), z.Not(), y]); M.AddBoolOr([y.Not(), z.Not(), x])
    for v in range(N):
        nb = sorted(linked[v])
        for i in range(len(nb)):
            for j in range(i + 1, len(nb)):
                add_triple(nb[i], v, nb[j])
    obj = []
    base = const
    for (a, b) in pairs:
        ee, nn = e.get((a, b), 0), ne[(a, b)]
        # cost = ee*(1-T) + nn*T
        base += ee
        obj.append((nn - ee) * t(a, b))
    M.Minimize(sum(obj))
    it = 0
    t0 = time.time()
    while True:
        it += 1
        S = cp_model.CpSolver()
        S.parameters.num_workers = workers
        S.parameters.max_time_in_seconds = tlimit
        st = S.Solve(M)
        if st not in (cp_model.OPTIMAL,):
            return None, st
        Tm = np.zeros((N, N), dtype=bool)
        for (a, b), var in T.items():
            if S.Value(var):
                Tm[a, b] = Tm[b, a] = True
        viol = []
        for v in range(N):
            nbr = np.nonzero(Tm[v])[0]
            if len(nbr) < 2: continue
            sub = Tm[np.ix_(nbr, nbr)]
            for (i, j) in np.argwhere(~sub):
                if i < j:
                    viol.append((nbr[i], v, nbr[j]))
        cost = base + int(round(S.ObjectiveValue()))
        if verbose: print('   it', it, 'relaxed cost', cost, 'viol', len(viol), 'triples', len(seen), round(time.time() - t0, 1), flush=True)
        if not viol:
            lab = [-1] * N; k = 0
            for v in range(N):
                if lab[v] < 0:
                    lab[v] = k
                    for w in np.nonzero(Tm[v])[0]: lab[w] = k
                    k += 1
            return cost, lab
        for (a, b, c) in viol[:50000]:
            add_triple(int(a), int(b), int(c))

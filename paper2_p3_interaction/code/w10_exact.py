"""EXACT decision  OPT(G) <= |H| + l  (all pairs may be edited) by lazy-triple SAT (CaDiCaL via pysat) with cardinality bound on the excess.
Variables T(u,v) = u,v in one cluster, for ALL pairs.  Excess tokens: hard-pair edits (pairs in no member) + (excess of each member: state O -> 1, X -> 2).
Triple (transitivity) clauses are added lazily (CEGAR): initial = all triples with two linked pairs (soft pair or edge); then violated triples of each model.
SAT answer is certified by checking the cluster graph and cost; UNSAT is a proof (relaxation is UNSAT)."""
import sys, os, itertools, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w8c_lib import *
from pysat.card import CardEnc, EncType


def exact_decide(n, E, members, l, verbose=True, max_iter=100000, add_per_iter=20000, timeout=None):
    t0 = time.time()
    soft = set()
    for t in members:
        soft |= pairs_of(t)
    pairs = list(itertools.combinations(range(n), 2))
    var = {}
    for i, pq in enumerate(pairs):
        var[pq] = i + 1
    nv = len(pairs)
    def T(u, v):
        return var[key(u, v)]
    s = Cadical153()
    nclauses = 0
    def add(c):
        nonlocal nclauses
        s.add_clause(c); nclauses += 1
    # linked adjacency (soft pairs and edges)
    adj = [set() for _ in range(n)]
    for pq in soft:
        adj[pq[0]].add(pq[1]); adj[pq[1]].add(pq[0])
    for pq in E:
        adj[pq[0]].add(pq[1]); adj[pq[1]].add(pq[0])
    seen = set()
    def add_triple(u, v, w):
        t = tuple(sorted((u, v, w)))
        if t in seen:
            return
        seen.add(t)
        a, b, c = T(t[0], t[1]), T(t[1], t[2]), T(t[0], t[2])
        add([-a, -b, c]); add([-a, -c, b]); add([-b, -c, a])
    for v in range(n):
        nb = sorted(adj[v])
        for i in range(len(nb)):
            for j in range(i + 1, len(nb)):
                add_triple(nb[i], v, nb[j])
    # excess tokens
    tokens = []
    for (a, b, c) in members:
        ab, bc, ac = T(a, b), T(b, c), T(a, c)
        q1 = nv + 1; q2 = nv + 2; nv += 2
        # O: not ab, not bc, not ac -> q1 ; X: not ab, not bc, ac -> q1,q2
        add([ab, bc, ac, q1]); add([ab, bc, -ac, q1]); add([ab, bc, -ac, q2])
        tokens += [q1, q2]
    for pq in pairs:
        if pq in soft:
            continue
        tokens.append(-T(*pq) if pq in E else T(*pq))
    enc = CardEnc.atmost(lits=tokens, bound=l, top_id=nv, encoding=EncType.seqcounter)
    for c in enc.clauses:
        add(c)
    if verbose:
        print('  exact: n', n, 'pairs', len(pairs), 'tokens', len(tokens), 'initial clauses', nclauses, 'triples', len(seen), round(time.time() - t0, 1), flush=True)
    it = 0
    while True:
        it += 1
        if timeout and time.time() - t0 > timeout:
            return 'TIMEOUT', None
        if not s.solve():
            if verbose: print('  exact: UNSAT after', it, 'iterations,', len(seen), 'triples', round(time.time() - t0, 1), flush=True)
            s.delete()
            return False, None
        model = s.get_model()
        ms = np.zeros(len(pairs) + 1, dtype=bool)
        for lit in model[:len(pairs)]:
            if lit > 0: ms[lit] = True
        Tm = np.zeros((n, n), dtype=bool)
        for (u, v), i in var.items():
            if ms[i]:
                Tm[u, v] = Tm[v, u] = True
        # violated triples: T[u,v] & T[v,w] & ~T[u,w]
        viol = []
        for v in range(n):
            nbr = np.nonzero(Tm[v])[0]
            if len(nbr) < 2: continue
            sub = Tm[np.ix_(nbr, nbr)]
            bad = np.argwhere(~sub)
            for (i, j) in bad:
                if i < j:
                    viol.append((nbr[i], v, nbr[j]))
                    if len(viol) >= add_per_iter: break
            if len(viol) >= add_per_iter: break
        if not viol:
            lab = labels_from(Tm)
            cost = partition_cost(n, E, lab)
            if verbose: print('  exact: SAT after', it, 'iterations; cost', cost, '|H| =', len(members), 'excess', cost - len(members), round(time.time() - t0, 1), flush=True)
            s.delete()
            return True, (lab, cost)
        for (u, v, w) in viol:
            add_triple(u, v, w)
        if verbose and it % 10 == 0:
            print('  exact: iter', it, 'triples', len(seen), round(time.time() - t0, 1), flush=True)


def labels_from(Tm):
    n = Tm.shape[0]
    lab = [-1] * n
    k = 0
    for v in range(n):
        if lab[v] < 0:
            lab[v] = k
            for w in np.nonzero(Tm[v])[0]:
                lab[w] = k
            k += 1
    return lab

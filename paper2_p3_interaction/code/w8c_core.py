"""Wave 8: SAT/CEGIS search for a K-core: a tight instance (OPT=|H|) with K clean anchors (ends of members,
no hard edge) such that every (K-1)-subset can be made singletons simultaneously by a ZE solution but all K cannot.

usage: python3 w8c_core.py n M K cmax [seconds] [seed]
Exact within the stated model: the search is complete for (n vertices, exactly M members, c<=cmax, anchors = vertices 0..K-1,
each anchor an end of exactly one member and incident to no hard edge)."""
import sys, os, time, random, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w8c_lib import *
from pysat.solvers import Cadical153
from pysat.card import CardEnc, EncType


def build(n, M, K, cmax, extra=None):
    vp = VarPool()
    P = list(itertools.combinations(range(n), 2))
    E = {p: vp.new() for p in P}
    Mem = {}
    for b in range(n):
        for a, c in itertools.combinations([v for v in range(n) if v != b], 2):
            if b < K:  # anchors are never centres
                continue
            Mem[(b, a, c)] = vp.new()
    cl = []
    cover = {p: [] for p in P}
    inmem = {v: [] for v in range(n)}
    for (b, a, c), m in Mem.items():
        cl += [[-m, E[key(a, b)]], [-m, E[key(b, c)]], [-m, -E[key(a, c)]]]
        for p in (key(a, b), key(b, c), key(a, c)):
            cover[p].append(m)
        for v in (a, b, c):
            inmem[v].append(m)
    S = {}
    for p in P:
        S[p] = vp.new()
        cl.append([-S[p]] + cover[p])
        for m in cover[p]:
            cl.append([-m, S[p]])
        # at most one cover
        for m1, m2 in itertools.combinations(cover[p], 2):
            cl.append([-m1, -m2])
    H0, H1 = {}, {}
    for p in P:
        H0[p] = vp.new(); H1[p] = vp.new()
        cl += [[-H0[p], -S[p]], [-H0[p], -E[p]], [H0[p], S[p], E[p]]]
        cl += [[-H1[p], -S[p]], [-H1[p], E[p]], [H1[p], S[p], -E[p]]]
    allm = list(Mem.values())
    enc = CardEnc.equals(lits=allm, bound=M, top_id=vp.n, encoding=EncType.seqcounter)
    vp.n = max(vp.n, enc.nv)
    cl += enc.clauses
    for v in range(n):
        enc = CardEnc.atmost(lits=inmem[v], bound=cmax, top_id=vp.n, encoding=EncType.seqcounter)
        vp.n = max(vp.n, enc.nv)
        cl += enc.clauses
    # anchors: in exactly one member, clean
    for x in range(K):
        enc = CardEnc.equals(lits=inmem[x], bound=1, top_id=vp.n, encoding=EncType.seqcounter)
        vp.n = max(vp.n, enc.nv)
        cl += enc.clauses
        for v in range(n):
            if v != x:
                cl.append([-E[key(x, v)], S[key(x, v)]])
    # symmetry breaking (light): anchors' centres increasing
    # feasibility scenarios
    scen = []
    for drop in range(K):
        sing = [x for x in range(K) if x != drop]
        c2, T = ze_cnf_cand(n, E, Mem, S, vp, sing)
        cl += c2
    return vp, E, Mem, S, H0, H1, cl


def ze_cnf_cand(n, E, Mem, S, vp, singletons):
    T = {}
    for p in itertools.combinations(range(n), 2):
        T[p] = vp.new()
    cl = []
    for u, v, w in itertools.combinations(range(n), 3):
        a, b, c = T[(u, v)], T[(v, w)], T[(u, w)]
        cl += [[-a, -b, c], [-a, -c, b], [-b, -c, a]]
    for p in T:
        # hard pair: T == E  when not soft
        cl += [[S[p], -T[p], E[p]], [S[p], T[p], -E[p]]]
    for (b, a, c), m in Mem.items():
        e1, e2, e3 = -T[key(a, b)], -T[key(b, c)], T[key(a, c)]
        cl += [[-m, e1, e2, e3], [-m, -e1, -e2], [-m, -e1, -e3], [-m, -e2, -e3]]
    for x in singletons:
        for v in range(n):
            if v != x:
                cl.append([-T[key(x, v)]])
    return cl, T


def decode(model_set, n, E, Mem):
    Eset = set(p for p, v in E.items() if v in model_set)
    members = [k_ for k_, v in Mem.items() if v in model_set]
    members = [(a, b, c) for (b, a, c) in members]
    return Eset, members


def run(n, M, K, cmax, seconds, seed):
    rng = random.Random(seed)
    vp, E, Mem, S, H0, H1, cl = build(n, M, K, cmax)
    solver = Cadical153(bootstrap_with=cl)
    t0 = time.time()
    it = 0
    nonanch = list(range(K, n))
    anch = list(range(K))
    while time.time() - t0 < seconds:
        if not solver.solve():
            print('UNSAT after', it, 'iterations: no K-core in this model', flush=True)
            return None
        model = set(l for l in solver.get_model() if l > 0)
        Eset, members = decode(model, n, E, Mem)
        lab = ze_solve(n, Eset, members, singletons=anch)
        it += 1
        if lab is None:
            print('FOUND candidate after', it, 'iterations', flush=True)
            return Eset, members
        # counterexample clause (+ permuted copies)
        for rep in range(40):
            if rep == 0:
                perm = list(range(n))
            else:
                a2 = anch[:]
                rng.shuffle(a2)
                o2 = nonanch[:]
                rng.shuffle(o2)
                perm = [0] * n
                for i, v in enumerate(anch):
                    perm[v] = a2[i]
                for i, v in enumerate(nonanch):
                    perm[v] = o2[i]
            lab2 = [None] * n
            for v in range(n):
                lab2[perm[v]] = lab[v]
            same = lambda u, v: lab2[u] == lab2[v]
            c = []
            for p in itertools.combinations(range(n), 2):
                if same(*p):
                    c.append(H0[p])
                else:
                    c.append(H1[p])
            for (b, a, cc), m in Mem.items():
                cnt = (0 if same(a, b) else 1) + (0 if same(b, cc) else 1) + (1 if same(a, cc) else 0)
                if cnt != 1:
                    c.append(m)
            solver.add_clause(c)
        if it % 200 == 0:
            print('iter', it, 'time', round(time.time() - t0, 1), flush=True)
    print('TIMEOUT after', it, 'iterations', flush=True)
    return 'timeout'


if __name__ == '__main__':
    n, M, K, cmax = map(int, sys.argv[1:5])
    sec = float(sys.argv[5]) if len(sys.argv) > 5 else 600
    seed = int(sys.argv[6]) if len(sys.argv) > 6 else 1
    res = run(n, M, K, cmax, sec, seed)
    if res and res != 'timeout':
        Eset, members = res
        print('n', n, 'edges', sorted(Eset))
        print('members', members)
        ok, why = check_instance(n, Eset, members)
        print('valid', ok, why)
        lab = ze_solve(n, Eset, members)
        print('tight ZE', lab is not None, lab)
        for drop in range(K):
            sing = [x for x in range(K) if x != drop]
            print('lock all but', drop, '->', ze_solve(n, Eset, members, singletons=sing) is not None)
        print('lock all ->', ze_solve(n, Eset, members, singletons=list(range(K))) is not None)

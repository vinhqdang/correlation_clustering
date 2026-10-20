"""Audit 8: Proposition 9 (tunnel-release model).  Tunnels = stubs (private pendant vertex) or paired hard K2 between two clean anchors.
f(R)=1 iff G minus tunnel edges of R has a zero-excess solution.  Checks: (i) monotone, (ii) f(R)=1 iff R hits all minimal non-lockable sets,
(iii) min |R| = TR optimum (brute force over all partitions of the tunnelled graph with only tunnel cuts as hard edits and one edit per member),
and the bounded search tree algorithm.  Also records whether unrestricted OPT - |H| can be smaller than the TR optimum (model caveat)."""
import sys, random, itertools
import numpy as np
from a8_lib import *
rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
N_INST = int(sys.argv[2]) if len(sys.argv) > 2 else 200
stats = {}; viol = []
def bump(k, v=1): stats[k] = stats.get(k, 0) + v

def clean_anchors(I):
    return [x for x in range(I.n) if not I.hard_nbrs(x) and any(x in t for t in I.H)]

def build(I, tunnels):
    """tunnels: list of ('stub', x) or ('pair', x1, x2). returns A of tunnelled graph and tunnel edge list"""
    ns = sum(1 for t in tunnels if t[0] == 'stub'); n2 = I.n + ns
    A = np.zeros((n2, n2), np.int8); A[:I.n, :I.n] = I.A; edges = []; nxt = I.n
    for t in tunnels:
        if t[0] == 'stub': A[t[1], nxt] = A[nxt, t[1]] = 1; edges.append((t[1], nxt)); nxt += 1
        else: A[t[1], t[2]] = A[t[2], t[1]] = 1; edges.append((t[1], t[2]))
    return A, edges

def f_of(A, edges, H, R, cache):
    key_ = frozenset(R)
    if key_ in cache: return cache[key_]
    A2 = A.copy()
    for i in R: u, v = edges[i]; A2[u, v] = A2[v, u] = 0
    J = Inst(A2, H); r = J.opt() == len(H); cache[key_] = r; return r

def run(I):
    anchors = clean_anchors(I)
    if len(anchors) < 2: return
    nt = rng.choice([2, 3]); tunnels = []; used = set(); pool = anchors[:]; rng.shuffle(pool)
    stubs_budget = 10 - I.n
    for _ in range(nt):
        if not pool: break
        x = pool.pop()
        if rng.random() < 0.5 or not pool:
            if stubs_budget <= sum(1 for t in tunnels if t[0] == 'stub'): continue
            tunnels.append(('stub', x))
        else:
            # find partner with disjoint Q, no shared member, nonadjacent
            cand = [y for y in pool if not I.A[x, y] and key(x, y) not in I.soft and not (I.Q(x) & I.Q(y))]
            if not cand: tunnels.append(('stub', x)) if stubs_budget > sum(1 for t in tunnels if t[0] == 'stub') else None; continue
            y = rng.choice(cand); pool.remove(y); tunnels.append(('pair', x, y))
    if len(tunnels) < 2: return
    A, edges = build(I, tunnels); J = Inst(A, I.H)
    if not J.valid(): viol.append('invalid'); return
    T = list(range(len(tunnels))); cache = {}
    f = {frozenset(R): f_of(A, edges, I.H, R, cache) for r in range(len(T) + 1) for R in itertools.combinations(T, r)}
    # (i) monotone
    for R in f:
        for R2 in f:
            if R <= R2 and f[R] and not f[R2]: viol.append(('mono', I.A.tolist(), I.H, tunnels)); return
    lock = {S: f[frozenset(T) - S] for S in f}          # lockable S
    cores = [S for S in f if not lock[S] and all(lock[S - {t}] for t in S)]
    for R in f:   # (ii)
        if f[R] != all(R & c for c in cores): viol.append(('hit', I.A.tolist(), I.H, tunnels)); return
    h = min([len(R) for R in f if f[R]], default=None)   # (iii) min hitting set
    # TR optimum by brute force over all partitions of J: hard edits only tunnel cuts, exactly one edit per member
    tun_pairs = {key(*e) for e in edges}
    mem_idx = [[J.pidx[p] for p in mpairs(t)] for t in J.H]
    ok = np.ones(len(J.cost), bool)
    for ix in mem_idx: ok &= (J.diff[:, ix].sum(1) == 1)
    for p in J.P:
        if p in J.soft or p in tun_pairs: continue
        ok &= ~J.diff[:, J.pidx[p]]
    cuts = sum(J.diff[:, J.pidx[p]].astype(int) for p in tun_pairs)
    tr = int(cuts[ok].min()) if ok.any() else None
    bump('inst'); 
    if h != tr: viol.append(('TRopt', I.A.tolist(), I.H, tunnels, h, tr)); return
    allopt = J.opt() - len(J.H)
    if h is not None and allopt < h: bump('unrestricted_OPT_below_TR_opt')
    if h is not None and allopt > h: viol.append('TR opt below true opt?!'); return
    bump('h=%s' % h); bump('maxcore%d' % max([len(c) for c in cores], default=0))
    if any(len(R) for R in f if not f[R]): bump('nonlockable_nontrivial')
    # bounded search tree with greedy core finding, calling f (the l=0 DP)
    def alg(R):
        R = set(R)
        if f_of(A, edges, I.H, R, cache): return len(R)
        S = set(T) - R
        for t in sorted(S):
            if not lock_of(S - {t}, R): S = S - {t}
        # S is now a minimal non-lockable subset of T minus R (within current R): branch
        best = None
        for t in S:
            r = alg(R | {t}); best = r if best is None or r < best else best
        return best
    def lock_of(S, R):   # lockable w.r.t. current R: remove all tunnels outside S
        return f_of(A, edges, I.H, set(T) - S, cache)
    a = alg(set())
    if a != h: viol.append(('alg', I.A.tolist(), I.H, tunnels, a, h)); return

n_ok = 0; tried = 0
while n_ok < N_INST and tried < 400 * N_INST:
    tried += 1
    n = rng.choice([6, 7, 7])
    A = gen_planted(n, rng, nedit=rng.randint(1, 4)) if rng.random() < 0.85 else gen_random(n, rng)
    H = greedy_pack(A, rng)
    if not H: continue
    I = Inst(A, H)
    if not I.valid(): continue
    b = stats.get('inst', 0); run(I); n_ok += stats.get('inst', 0) - b
    if viol: break
print(stats); print('VIOLATIONS:', viol[:2])

"""Audit 8: Lemma 8 (EX(a) fan-out, NF(gamma) frozen components, bounded-degree Grid Tiling algorithm)."""
import sys, random, itertools
import numpy as np
from a8_lib import *
rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
N_INST = int(sys.argv[2]) if len(sys.argv) > 2 else 300
stats = {}; viol = []
def bump(k, v=1): stats[k] = stats.get(k, 0) + v

def hard_comps(I):
    n = I.n; par = list(range(n))
    def f(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    for u, v in itertools.combinations(range(n), 2):
        if I.A[u, v] and (u, v) not in I.soft: par[f(u)] = f(v)
    g = {}
    for v in range(n): g.setdefault(f(v), []).append(v)
    return list(g.values())

def test_inst(I):
    D = max(len(I.Q(v)) for v in range(I.n)); nH = len(I.H)
    dh = [len(I.hard_nbrs(v)) for v in range(I.n)]
    comps = hard_comps(I)
    for l in (0, 1, 2):
        sel = np.nonzero(I.cost <= nH + l)[0]
        for k in sel:
            Upairs = [p for p in I.P if p not in I.soft and I.diff[k, I.pidx[p]]]
            u = len(Upairs)
            for (x, y) in Upairs:
                bump('EXa_checks')
                if max(dh[x], dh[y]) > D + u or dh[x] > D + l or dh[y] > D + l: viol.append(('EXa', I.A.tolist(), I.H, k)); return
        # NF(gamma): for OPTIMAL partitions with cost <= nH + l
        if I.opt() <= nH + l:
            B = 2 * D + 3 * l; opt = I.opt()
            for k in np.nonzero(I.cost == opt)[0]:
                for comp in comps:
                    if max(dh[v] for v in comp) > B:
                        bump('NFg_checks')
                        cl = {v for v in range(I.n) if I.lab[k, v] == I.lab[k, comp[0]]}
                        if cl != set(comp): viol.append(('NFg', I.A.tolist(), I.H, k)); return
                        # frozen = no HARD pair of it edited
                        if any(I.diff[k, I.pidx[key(u, v)]] for u, v in itertools.combinations(comp, 2) if key(u, v) not in I.soft): viol.append(('NFg-hard',)); return
                        # paper wording: "no pair of it is edited" (soft pairs too)
                        if any(I.diff[k, I.pidx[key(u, v)]] for u, v in itertools.combinations(comp, 2)): bump('NFg_soft_pair_edited_inside_frozen')

n_ok = 0; tried = 0
while n_ok < N_INST and tried < 100 * N_INST:
    tried += 1
    n = rng.choice([7, 8, 9, 10])
    r = rng.random()
    if r < 0.5: A = gen_planted(n, rng, nedit=rng.randint(1, 3), bigcluster=True)
    elif r < 0.8: A = gen_planted(n, rng)
    else: A = gen_random(n, rng)
    H = greedy_pack(A, rng)
    if not H: continue
    I = Inst(A, H)
    if not I.valid() or not I.is_pair_maximal(): continue
    n_ok += 1; bump('inst')
    test_inst(I)
    if viol: break
print(stats); print('VIOLATIONS:', viol[:2])

# explicit tiny example of the wording issue: K_9 minus edge ac, H = {(a,b,c)}
n = 9; A = np.ones((n, n), np.int8) - np.eye(n, dtype=np.int8); A[0, 2] = A[2, 0] = 0
I = Inst(A, [(0, 1, 2)]); D = max(len(I.Q(v)) for v in range(n))
ks = np.nonzero(I.cost == I.opt())[0]
print('example K9-ac: valid', I.valid(), 'pairmax', I.is_pair_maximal(), 'OPT', I.opt(), 'D', D, 'dh(3)=', len(I.hard_nbrs(3)),
      'B(l=0)=', 2 * D, 'optima', len(ks), 'soft pair ac inserted in optimum:', bool(I.diff[ks[0], I.pidx[(0, 2)]]))

# ---------- (b) bounded-degree Grid Tiling
def gt_brute(k, n, S):
    for a in itertools.product(range(n), repeat=k):
        for b in itertools.product(range(n), repeat=k):
            if all((a[i], b[j]) in S[i][j] for i in range(k) for j in range(k)): return True
    return False
def gt_alg(k, n, S, d):
    for a1 in range(n):
        bch = []
        for j in range(k):
            ch = [b for b in range(n) if (a1, b) in S[0][j]]
            if len(ch) > d: raise Exception('degree')
            bch.append(ch)
        for bs in itertools.product(*bch):
            ach = []
            for i in range(k):
                ch = [a for a in range(n) if (a, bs[0]) in S[i][0]]
                if len(ch) > d: raise Exception('degree')
                ach.append(ch)
            for as_ in itertools.product(*ach):
                if all((as_[i], bs[j]) in S[i][j] for i in range(k) for j in range(k)): return True
    return False
cnt = 0; yes = 0
for _ in range(3000):
    k = rng.choice([2, 3]); n = rng.choice([3, 4, 5]); d = rng.choice([1, 2])
    S = [[set() for _ in range(k)] for _ in range(k)]
    for i in range(k):
        for j in range(k):
            for _m in range(d):
                perm = list(range(n)); rng.shuffle(perm)
                for a_ in range(n):
                    if rng.random() < 0.85: S[i][j].add((a_, perm[a_]))
    ok = all(sum(1 for a_ in range(n) if (a_, b_) in S[i][j]) <= d and sum(1 for b2 in range(n) if (b_, b2) in S[i][j]) <= d
             for i in range(k) for j in range(k) for b_ in range(n))
    if not ok: continue
    cnt += 1; r1 = gt_brute(k, n, S); r2 = gt_alg(k, n, S, d); yes += r1
    assert r1 == r2
print('GridTiling: instances', cnt, 'yes', yes, 'alg==brute on all')

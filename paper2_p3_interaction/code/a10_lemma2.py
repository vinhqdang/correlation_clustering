"""Lemma 2: OPT_Uonly - |H| = min{|U| : Sigma(G xor U) nonempty}  (identity, independent brute force), plus OPT vs U-only; plus locality claim."""
import sys, random, itertools
import numpy as np
from a10_lib import *

def rand_inst(n, rng, mmax):
    while True:
        p = rng.uniform(0.3, 0.7)
        A = np.zeros((n, n), bool)
        for i in range(n):
            for j in range(i + 1, n):
                if rng.random() < p: A[i, j] = A[j, i] = True
        P3 = [(a, b, c) for b in range(n) for a in range(n) for c in range(a + 1, n)
              if len({a, b, c}) == 3 and A[a, b] and A[b, c] and not A[a, c]]
        rng.shuffle(P3); H = []; used = set()
        for t in P3:
            ps = pairs_of(t)
            if any(x in used for x in ps): continue
            used.update(ps); H.append(t)
            if len(H) >= mmax: break
        if H: return A, H

def uonly_opt(A, H):
    n = len(A); P = partitions(n)
    same = P[:, :, None] == P[:, None, :]; diff = same != A[None]
    iu = np.triu_indices(n, 1); D = diff[:, iu[0], iu[1]]
    idx = {(int(i), int(j)): k for k, (i, j) in enumerate(zip(*iu))}
    pk = lambda x, y: idx[(min(x, y), max(x, y))]
    cost = D.sum(1); ok = np.ones(len(P), bool)
    for (a, b, c) in H:
        ok &= D[:, [pk(a, b), pk(b, c), pk(a, c)]].sum(1) == 1
    return int(cost[ok].min()) if ok.any() else None, int(cost.min())

def min_U(A, H, maxu=3):
    n = len(A); o = owner_map(H)
    hard = [(i, j) for i in range(n) for j in range(i + 1, n) if frozenset((i, j)) not in o]
    for u in range(maxu + 1):
        for U in itertools.combinations(hard, u):
            B = A.copy()
            for x, y in U: B[x, y] = B[y, x] = not B[x, y]
            sg, _ = sigma_bruteforce(B, H)
            if sg: return u
    return None

rng = random.Random(int(sys.argv[1])); N = int(sys.argv[2])
st = dict(tot=0, ident_ok=0, ident_bad=0, opt_lt_uonly=0, opt_eq=0, none=0)
for it in range(N):
    n = rng.randint(5, 7); A, H = rand_inst(n, rng, rng.randint(1, 3))
    uo, opt = uonly_opt(A, H)
    if uo is not None and uo-len(H)>4: continue
    mu = min_U(A, H, (uo-len(H)) if uo is not None else 3)
    st['tot'] += 1
    if (uo is None) != (mu is None):
        st['xor'] = st.get('xor',0)+1; print('XOR', uo, mu, H)
    if uo is None or mu is None: st['none'] += 1; continue
    if uo - len(H) == mu: st['ident_ok'] += 1
    else: st['ident_bad'] += 1; print("IDENT FAIL", A.astype(int).tolist(), H, uo, mu)
    if opt < uo: st['opt_lt_uonly'] += 1; ex = (A.copy(), H, opt, uo)
    else: st['opt_eq'] += 1
print(st)
if st['opt_lt_uonly']:
    A, H, opt, uo = ex; print("example OPT<Uonly:", [(i, j) for i in range(len(A)) for j in range(i+1, len(A)) if A[i, j]], H, "OPT", opt, "Uonly", uo)

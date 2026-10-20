"""Flip-gadget search with clean private sockets: members I1, I2 each have a private vertex (in no other member, no hard edge) s1, s2.
Then test the chain composition N=2 (identify s2 of gadget 0 with s1 of gadget 1) by ILP: tight, and relation flip.
n = 9 (anchor = vertex 8)."""
import sys, itertools, random, numpy as np, time, pickle
from w6c_lib import *
from w6c_compose2 import sfix, rel
n = int(sys.argv[1]); m = int(sys.argv[2]); T = float(sys.argv[3]); seed = int(sys.argv[4])
rng = random.Random(seed)
lab, same, P = part_table(n); same_f = same.astype(np.float32); pidx = {p: i for i, p in enumerate(P)}
def sig(mem):
    a, b, c = mem
    s = same[:, [pidx[key(a, b)], pidx[key(b, c)], pidx[key(a, c)]]].astype(int)
    out = np.full(len(s), -1)
    for code, pat in enumerate(((0, 1, 0), (1, 0, 0), (1, 1, 1))): out[(s == np.array(pat)).all(1)] = code
    return out
def rand_config():
    for _ in range(300):
        x = n - 1; mems = []
        y, z = rng.sample(range(n - 1), 2)
        P_ = rng.choice([[x, y, z], [x, z, y], [y, z, x], [z, y, x]])
        mems.append(tuple(P_))
        for _ in range(m - 1):
            for _ in range(60):
                t = tuple(rng.sample(range(n - 1), 3))
                if all(len(set(t) & set(u)) <= 1 for u in mems): mems.append(t); break
            else: break
        if len(mems) == m: return mems
    return None
t0 = time.time(); stats = dict(cfg=0, flip=0, comp=0); found = []
while time.time() - t0 < T:
    mems = rand_config()
    if mems is None: continue
    cnt = {}
    for mm in mems:
        for v in mm: cnt[v] = cnt.get(v, 0) + 1
    priv = {i: [v for v in mems[i] if cnt[v] == 1] for i in range(1, m)}
    if sum(1 for i in priv if priv[i]) < 2: continue
    stats['cfg'] += 1
    soft = set(); fe = set()
    for mm in mems: soft |= pairs_of(mm); fe |= {key(mm[0], mm[1]), key(mm[1], mm[2])}
    x = n - 1
    hard = [p for p in P if p not in soft and x not in p]
    h = len(hard)
    if h > 14: continue
    S = np.stack([sig(mm) for mm in mems], 1); valid = (S >= 0).all(1)
    sing = np.ones(len(same), bool)
    for t in range(n):
        if t != x: sing &= ~same[:, pidx[key(x, t)]]
    CH = 1 << min(h, 11)
    for c0 in range(0, 1 << h, CH):
        js = np.arange(c0, min(1 << h, c0 + CH))
        Hm = ((js[:, None] >> np.arange(h)[None, :]) & 1).astype(np.float32)
        E = np.zeros((len(js), len(P)), np.float32)
        for q in fe: E[:, pidx[q]] = 1
        for i, q in enumerate(hard): E[:, pidx[q]] = Hm[:, i]
        cost = same_f @ (1 - 2 * E.T) + E.sum(1)[None, :]
        okm = (cost == m)
        for jj in range(len(js)):
            ok = okm[:, jj] & valid
            if not ok.any(): continue
            j = int(js[jj])
            hv = {v for i, q in enumerate(hard) if (j >> i) & 1 for v in q}
            for i1, i2 in itertools.permutations(range(1, m), 2):
                s1c = [v for v in priv[i1] if v not in hv]; s2c = [v for v in priv[i2] if v not in hv]
                if not s1c or not s2c: continue
                Rl = {(s[i1], s[i2]) for s in S[ok & sing]}
                Rr = {(s[i1], s[i2]) for s in S[ok & ~sing]}
                if len(Rl) == 2 and len({a for a, _ in Rl}) == 2 and len({b for _, b in Rl}) == 2:
                    D1 = {a for a, _ in Rl}; D2 = {b for _, b in Rl}; f = dict(Rl)
                    new = [(d, e) for d, e in Rr if d in D1 and e in D2 and f[d] != e]
                    if not new: continue
                    stats['flip'] += 1
                    rec = (mems, j, i1, i2, [tuple(map(int, t)) for t in sorted(Rl)], [tuple(map(int, t)) for t in sorted(Rr)], [tuple(map(int, t)) for t in new], hard)
                    from w6c_compose3 import build as cbuild
                    for s1 in s1c:
                        for s2 in s2c:
                            A, H, xs, I1f, I2l = cbuild(rec, 2, s1, s2)
                            ok_, why = check_packing(A, H)
                            if not ok_: continue
                            if ilp_opt(A)[0] != len(H): continue
                            stats['comp'] += 1
                            allL = rel(A, H, xs, xs, I1f, I2l)
                            rels = [rel(A, H, xs, [x_ for x_ in xs if x_ != xs[g]], I1f, I2l) for g in range(2)]
                            flipc = len(allL) >= 1 and all(any(p not in allL for p in R) for R in rels)
                            print('COMPOSITE tight', mems, j, (i1, i2), (s1, s2), 'allL', sorted(allL), 'rels', [sorted(R) for R in rels], 'FLIPCHAIN' if flipc else '', flush=True)
                            found.append((rec, s1, s2, flipc))
                            pickle.dump(found, open(outpath(f'flip5_{n}_{m}_{seed}.pkl'), 'wb'))
print('n', n, 'm', m, stats, flush=True)

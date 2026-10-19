"""Search gadgets with an anchor x (member P, lock) and two clean socket vertices l, r (ends of members, no hard edge at l, r).
flag f(v)=1 iff v is a singleton cluster.  Spec (GR): zero-excess solutions with x singleton never have f(l)=f(r)=1, some have f(l)=1,f(r)=0
(and some f(l)=0,f(r)=1);  some zero-excess solution with x non-singleton has f(l)=f(r)=1."""
import sys, itertools, random, numpy as np, time, pickle
from w6c_lib import *
n = int(sys.argv[1]); m = int(sys.argv[2]); T = float(sys.argv[3]); seed = int(sys.argv[4])
rng = random.Random(seed)
lab, same, P = part_table(n); same_f = same.astype(np.float32); pidx = {p: i for i, p in enumerate(P)}
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
t0 = time.time(); found = []; stats = dict(cfg=0, pat=0)
while time.time() - t0 < T:
    mems = rand_config()
    if mems is None: continue
    stats['cfg'] += 1
    soft = set(); fe = set(); inm = {}
    for i, mm in enumerate(mems):
        soft |= pairs_of(mm); fe |= {key(mm[0], mm[1]), key(mm[1], mm[2])}
        for pos, v in enumerate(mm): inm.setdefault(v, []).append((i, pos))
    x = n - 1
    # socket candidates: vertices (not x) that are an END in every member containing them
    cands = [v for v in range(n - 1) if v in inm and all(pos in (0, 2) for (_, pos) in inm[v])]
    if len(cands) < 2: continue
    hard = [p for p in P if p not in soft and x not in p]
    h = len(hard)
    if h > 13: continue
    Hm = ((np.arange(2 ** h)[:, None] >> np.arange(h)[None, :]) & 1).astype(np.float32)
    E = np.zeros((2 ** h, len(P)), np.float32)
    for p in fe: E[:, pidx[p]] = 1
    for i, p in enumerate(hard): E[:, pidx[p]] = Hm[:, i]
    cost = same_f @ (1 - 2 * E.T) + E.sum(1)[None, :]
    okm = cost == m
    # singleton indicator per vertex
    sing = {}
    for v in range(n):
        s = np.ones(len(same), bool)
        for t in range(n):
            if t != v: s &= ~same[:, pidx[key(v, t)]]
        sing[v] = s
    for j in range(2 ** h):
        ok = okm[:, j]
        if not ok.any(): continue
        stats['pat'] += 1
        hv = {v for i, p in enumerate(hard) if (j >> i) & 1 for v in p}
        L = ok & sing[x]; Rr = ok & ~sing[x]
        if not L.any() or not Rr.any(): continue
        for l, r in itertools.permutations([c for c in cands if c not in hv], 2):
            if l > r: continue
            fl, fr = sing[l], sing[r]
            if (L & fl & fr).any(): continue
            if not (L & fl & ~fr).any() or not (L & ~fl & fr).any(): continue
            if not (Rr & fl & fr).any(): continue
            found.append((mems, j, l, r, hard)); print('SOCKET-GADGET', mems, j, 'sockets', l, r, flush=True)
            pickle.dump(found, open(f'/tmp/claude-0/scratch/sock_{n}_{m}_{seed}.pkl', 'wb'))
print('n', n, 'm', m, stats, 'found', len(found), flush=True)

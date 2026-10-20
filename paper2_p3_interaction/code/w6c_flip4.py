"""Random search for a persistent-bit 'flip' lock gadget with an arbitrary number of members (default 4) on n<=8 vertices.
Anchor member P=(x,y,z) (x clean private end).  Interface members I1,I2 (any two other members).
R_lock = projection on (s_I1,s_I2) of zero-excess solutions with x singleton; R_rel = same for x non-singleton.
FLIP: R_lock is a bijection between two 2-sets D1->D2 and R_rel contains a pair (d,e), d in D1, e in D2, e != R_lock(d)."""
import sys, itertools, random, numpy as np, time
from w6c_lib import *
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
    for _ in range(200):
        verts = list(range(n)); x = n - 1
        mems = []
        # anchor member: x at end, two others random
        y, z = rng.sample(range(n - 1), 2)
        P_ = [x, y, z] if rng.random() < .5 else [z, y, x]
        if rng.random() < .5: P_ = [x, z, y] if P_[0] == x else [y, z, x]
        mems.append(tuple(P_))
        ok = True
        for _ in range(m - 1):
            for _ in range(50):
                t = tuple(rng.sample(range(n - 1), 3))
                if all(len(set(t) & set(u)) <= 1 for u in mems): mems.append(t); break
            else: ok = False; break
        if ok and len(mems) == m: return mems
    return None
t0 = time.time(); tested = 0; found = []
stats = {'cfg': 0, 'bij': 0}
while time.time() - t0 < T:
    mems = rand_config()
    if mems is None: continue
    stats['cfg'] += 1
    soft = set(); fe = set()
    for mm in mems: soft |= pairs_of(mm); fe |= {key(mm[0], mm[1]), key(mm[1], mm[2])}
    x = mems[0][0] if mems[0][0] == n - 1 else mems[0][2]
    hard = [p for p in P if p not in soft and x not in p]
    h = len(hard)
    if h > 13: continue
    Hm = ((np.arange(2 ** h)[:, None] >> np.arange(h)[None, :]) & 1).astype(np.float32)
    E = np.zeros((2 ** h, len(P)), np.float32)
    for p in fe: E[:, pidx[p]] = 1
    for i, p in enumerate(hard): E[:, pidx[p]] = Hm[:, i]
    cost = same_f @ (1 - 2 * E.T) + E.sum(1)[None, :]
    S = np.stack([sig(mm) for mm in mems], 1); valid = (S >= 0).all(1)
    sing = np.ones(len(same), bool)
    for t in range(n):
        if t != x: sing &= ~same[:, pidx[key(x, t)]]
    for j in range(2 ** h):
        ok = (cost[:, j] == m) & valid
        if not ok.any(): continue
        tested += 1
        for i1, i2 in itertools.permutations(range(1, m), 2):
            Rl = {(s[i1], s[i2]) for s in S[ok & sing]}
            Rr = {(s[i1], s[i2]) for s in S[ok & ~sing]}
            if len(Rl) == 2 and len({a for a, _ in Rl}) == 2 and len({b for _, b in Rl}) == 2:
                stats['bij'] += 1
                D1 = {a for a, _ in Rl}; D2 = {b for _, b in Rl}; f = dict(Rl)
                new = [(d, e) for d, e in Rr if d in D1 and e in D2 and f[d] != e]
                if new:
                    found.append((mems, j, i1, i2, [tuple(map(int,t)) for t in sorted(Rl)], [tuple(map(int,t)) for t in sorted(Rr)], [tuple(map(int,t)) for t in new], hard)); import pickle; pickle.dump(found, open(outpath(f'flip_{n}_{m}_{seed}.pkl'),'wb'))
print('n', n, 'm', m, 'configs', stats['cfg'], 'patterns with solutions', tested, 'bijective locked relations', stats['bij'], 'FLIPs', len(found), flush=True)

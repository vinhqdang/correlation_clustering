"""Template-seeded socket-gadget search (n=9).  Vertices: x=8 anchor; template members P=(x,y,z), Q=(y,z',w), Ml=(l,p,z), Mr=(r,p',z')
plus k random extra members.  Spec as in w6c_sock.py (GR)."""
import sys, itertools, random, numpy as np, time, pickle
from w6c_lib import *
n = 9; k_extra = int(sys.argv[1]); T = float(sys.argv[2]); seed = int(sys.argv[3])
rng = random.Random(seed)
lab, same, P = part_table(n); same_f = same.astype(np.float32); pidx = {p: i for i, p in enumerate(P)}
y, z, zp, p_, pp, l, r, w, x = range(9)
base = [(x, y, z), (y, zp, w), (l, p_, z), (r, pp, zp)]
def run(mems, sockets):
    soft = set(); fe = set()
    for mm in mems: soft |= pairs_of(mm); fe |= {key(mm[0], mm[1]), key(mm[1], mm[2])}
    cl = set(sockets) | {x}
    hard = [q for q in P if q not in soft and not (set(q) & cl)]
    h = len(hard)
    if h > 15: return []
    out = []
    sing = {}
    for v in cl:
        s = np.ones(len(same), bool)
        for t in range(n):
            if t != v: s &= ~same[:, pidx[key(v, t)]]
        sing[v] = s
    CH = 1 << min(h, 11)
    for c0 in range(0, 1 << h, CH):
        js = np.arange(c0, min(1 << h, c0 + CH))
        Hm = ((js[:, None] >> np.arange(h)[None, :]) & 1).astype(np.float32)
        E = np.zeros((len(js), len(P)), np.float32)
        for q in fe: E[:, pidx[q]] = 1
        for i, q in enumerate(hard): E[:, pidx[q]] = Hm[:, i]
        cost = same_f @ (1 - 2 * E.T) + E.sum(1)[None, :]
        okm = cost == len(mems)
        for jj in range(len(js)):
            ok = okm[:, jj]
            if not ok.any(): continue
            L = ok & sing[x]; Rr = ok & ~sing[x]
            if not L.any() or not Rr.any(): continue
            fl, fr = sing[sockets[0]], sing[sockets[1]]
            if (L & fl & fr).any(): continue
            if not (L & fl & ~fr).any() or not (L & ~fl & fr).any(): continue
            if not (Rr & fl & fr).any(): continue
            out.append(int(js[jj]))
    return out, hard
t0 = time.time(); found = []; tested = 0
while time.time() - t0 < T:
    mems = list(base)
    ok = True
    for _ in range(k_extra):
        for _ in range(80):
            t = tuple(rng.sample(range(8), 3))
            if all(len(set(t) & set(u)) <= 1 for u in mems): mems.append(t); break
        else: ok = False
    if not ok: continue
    tested += 1
    res = run(mems, (l, r))
    if res and res[0]:
        found.append((mems, res[0][:3], res[1])); print('FOUND', mems, res[0][:3], flush=True)
        pickle.dump([(m_, j_[0], l, r, h_) for m_, j_, h_ in found], open(outpath(f'sock2_{k_extra}_{seed}.pkl'), 'wb'))
print('k_extra', k_extra, 'tested', tested, 'found', len(found), flush=True)

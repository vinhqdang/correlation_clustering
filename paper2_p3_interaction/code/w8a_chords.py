"""Path 0..2m (or 0..N) plus one or two chords: for every packing H of size >= nu (all are enumerated) compute far_count(H) = min over optima of the
number of hard edits with both ends covered and I(H)-distance >= 2.  Report min over H (best repack), max over H (worst packing), and #exchanges needed
from the worst packing to a best one."""
import sys, itertools, time, collections
from w8a_lib import *
from w8a_kappa import far_count
def scan(n, extra, maxpk=100000):
    E = [(i, i + 1) for i in range(n - 1)] + list(extra)
    g = G8(n, E)
    opt, sols = g.opt(limit=5000)
    optima = [edits_of_partition(g, part_from_mask(s, n)) for s in sols]
    nu, _ = nu_mod(g)
    packs = enum_packings(g, minsize=nu, maxcount=maxpk)
    vals = []
    for H in packs:
        H = list(H)
        if not pair_maximal(g, H): continue
        vals.append((far_count(H, optima), H))
    mn = min(v for v, _ in vals); mx = max(v for v, _ in vals)
    good = [H for v, H in vals if v == mn]
    # min exchanges from worst to a good one
    worst = [H for v, H in vals if v == mx]
    ex = min(sum(1 for t in W if t not in set(Hg)) for W in worst for Hg in good) if mx > mn else 0
    return opt, nu, len(vals), mn, mx, ex
if __name__ == '__main__':
    t0 = time.time()
    for n in (9, 11, 13):
        res = collections.Counter(); rows = []
        for i in range(n):
            for j in range(i + 3, n):
                opt, nu, npk, mn, mx, ex = scan(n, [(i, j)])
                res[(opt - nu, mn, mx, ex)] += 1
                if mx > 0 or mn > 0: rows.append(((i, j), opt, nu, npk, mn, mx, ex))
        print('path n=%d + one chord:  (OPT-nu, min far, max far, exchanges worst->good) counts' % n, dict(res), '%.0fs' % (time.time() - t0), flush=True)
        for r in rows[:6]: print('   ', r)

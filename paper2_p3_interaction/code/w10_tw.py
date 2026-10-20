"""Pathwidth upper bound of I(H) for the reduction instance by the explicit left-to-right sweep (interval supergraph):
members = intervals [min coord, max coord]; two members sharing a vertex intersect; pw(I(H)) <= max point cover - 1.
Tunnels are hard pairs and do not appear in I(H)."""
import sys, random, time
from w10_red import *

def coord(name):
    if name[0] in ('z', 'y', 'x') and isinstance(name[1], tuple) and name[1][0] in ('A', 'B'):
        return 2 * name[2] + 1
    if name[0] == 'c':
        return None  # computed from alpha below
    # base vertex (tag, thin, i) or (tag, thin) in thin mode
    thin_ = name[1]
    kind, idx = thin_
    return idx // 2 if kind == 's' else idx

def sweep_width(W, meta):
    n = meta['n']
    names = W.B.names
    sep = meta['sep_group']
    cd = {}
    for v in names:
        if v[0] == 'c':
            cd[v] = 2 * sep[v[3]] + 1
        else:
            cd[v] = coord(v)
    events = {}
    ivs = []
    for (a, b, c) in W.B.members:
        cs = [cd[names[x]] for x in (a, b, c)]
        ivs.append((min(cs), max(cs)))
    mx = max(h for l, h in ivs)
    cover = [0] * (mx + 2)
    for l, h in ivs:
        for t in range(l, h + 1): cover[t] += 1
    return max(cover) - 1, ivs

if __name__ == '__main__':
    for (k, n) in [(1, 2), (1, 3), (2, 2), (2, 3), (2, 4), (3, 2)]:
        rng = random.Random(1)
        S = [[set((a, b) for a in range(1, n + 1) for b in range(1, n + 1) if rng.random() < 0.6) or {(1, 1)} for _ in range(k)] for _ in range(k)]
        for p in (8,):
            W, meta = build_gt(k, n, S, p, thin_mode=False)
            I = Inst(W)
            pw, ivs = sweep_width(W, meta)
            print('k', k, 'n', n, 'p', p, 'vertices', I.n, 'members', len(I.members), 'c', I.c, 'm', meta['m'], 'pw(I(H)) <=', pw, flush=True)

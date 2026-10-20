"""independent check of the path decomposition of I(H): bags = members whose coordinate interval contains t; verify validity by brute force
(every pair of members sharing a vertex has intersecting intervals) and report max bag, max span, for several n."""
import sys, random
from a15_lib import *
from w10_tw import coord
def run(k, n, p):
    rng = random.Random(5)
    S = [[set((a, b) for a in range(1, n + 1) for b in range(1, n + 1) if rng.random() < 0.6) or {(1, 1)} for _ in range(k)] for _ in range(k)]
    W, meta = build_gt(k, n, S, p, thin_mode=False)
    I = Inst(W)
    names = W.B.names; sep = meta['sep_group']
    def cd(nm):
        if nm[0] == 'c': return 2 * sep[nm[3]] + 1
        return coord(nm)
    cv = [cd(nm) for nm in names]
    ivs = []
    for (a, b, c) in I.members:
        xs = (cv[a], cv[b], cv[c]); ivs.append((min(xs), max(xs)))
    span = max(h - l for l, h in ivs)
    # validity: members sharing a vertex
    byv = {}
    for mi, t in enumerate(I.members):
        for v in t: byv.setdefault(v, []).append(mi)
    bad = 0
    for v, ms in byv.items():
        lo = max(ivs[x][0] for x in ms); hi = min(ivs[x][1] for x in ms)
        if lo > hi: bad += 1       # Helly for intervals: pairwise intersect iff common point
    mx = max(h for l, h in ivs)
    cover = [0] * (mx + 2)
    for l, h in ivs:
        for t in range(l, h + 1): cover[t] += 1
    print('k', k, 'n', n, 'p', p, '|V|', I.n, '|H|', len(I.members), 'm', meta['m'], 'max span', span, 'bad vertices', bad, 'max bag', max(cover), '-> pw <=', max(cover) - 1, 'c', I.c, flush=True)
if __name__ == '__main__':
    p = int(sys.argv[1])
    for (k, n) in [(1, 2), (1, 3), (1, 4), (1, 5)]:
        run(k, n, p)

"""u_max = max over super-vertices of (#pad matchings + #host matchings) computed in the thin-guest model (no thick build)."""
import sys, random
from collections import Counter
from w10_red import *
for (k, n) in [(1, 2), (2, 2), (2, 3), (3, 2), (3, 3), (3, 4), (4, 2), (4, 3), (5, 2)]:
    rng = random.Random(1)
    S = [[set((a, b) for a in range(1, n + 1) for b in range(1, n + 1)) for _ in range(k)] for _ in range(k)]   # all tiles present: worst case
    W, meta = build_gt(k, n, S, 2, thin_mode=True)
    cnt = Counter()
    for f in W.free:
        for x in f:
            if isinstance(x, tuple) and len(x) == 2 and isinstance(x[1], tuple) and x[1][0] in ('s', 'e'):
                cnt[x] += 1
    # host matchings: each perm guest z hosts at ('s',4u+2): one matching per group per chain
    for tag, ch in meta['chains'].items():
        for u in range(meta['U']):
            cnt[(tag, ('s', 4 * u + 2))] += 1
    print('k', k, 'n', n, 'u_max', max(cnt.values()), 'argmax', max(cnt, key=cnt.get)[1], flush=True)

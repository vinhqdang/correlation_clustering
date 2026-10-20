"""Tree-like packings (each new member shares exactly one vertex with the earlier ones, so no soft triangle; hubs allowed, c up to 4),
random hard pairs. n<=9. Sigma by partition enumeration; closure under all 729 rules, 2-decomposability, own path consistency vs truth."""
import sys, random, itertools
import numpy as np
from a10_lib import *
from a10_closure import pc_feasible
rng = random.Random(int(sys.argv[1])); N = int(sys.argv[2])
st = dict(tot=0, nonempty=0, ge3=0, closed=0, dec2=0, pc_ok=0, pc_bad=0, c=dict())
for it in range(N):
    m = rng.randint(2, 4); n = 3 + 2 * (m - 1)
    H = []; nxt = 3
    first = [0, 1, 2]; rng.shuffle(first); H.append(tuple(first))
    hub = rng.random() < .4; hv = H[0][rng.randrange(3)]
    for _ in range(m - 1):
        old = hv if hub else rng.choice([v for t in H for v in t])
        tri = [old, nxt, nxt + 1]; nxt += 2; rng.shuffle(tri); H.append(tuple(tri))
    A = np.zeros((n, n), bool)
    for (a, b, c) in H: A[a, b] = A[b, a] = A[b, c] = A[c, b] = True
    o = owner_map(H)
    ok = True
    for (a, b, c) in H:
        # a-c must be non-edge: it is, unless another member made it an edge (impossible, disjoint pairs)
        pass
    pbias = rng.choice([.2, .5, .8])
    for x, y in itertools.combinations(range(n), 2):
        if frozenset((x, y)) not in o: A[x, y] = A[y, x] = rng.random() < pbias
    if not valid(A, H) or soft_triangles(n, H): continue
    sig, _ = sigma_bruteforce(A, H)
    st['tot'] += 1; c = cmax(n, H); st['c'][c] = st['c'].get(c, 0) + 1
    st['nonempty'] += bool(sig); st['ge3'] += len(sig) >= 3
    if closed_all_rules(sig): st['closed'] += 1
    else: print("NOT CLOSED", H, sorted(sig))
    st['dec2'] += two_decomposable(sig, len(H))
    if pc_feasible(A, H) == bool(sig): st['pc_ok'] += 1
    else: st['pc_bad'] += 1
print(st)

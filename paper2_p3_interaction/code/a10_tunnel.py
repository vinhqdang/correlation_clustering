"""Cor A3 / Thm C(i): in soft-triangle-free instances, with stubs at anchors, G+stubs(S) has a ZE solution iff Sigma!=0 and every pair (and single) of S is a face.
Direct brute force on G+stubs (all partitions) vs prediction from faces of G."""
import sys, random, itertools
import numpy as np
from a10_lib import *
rng = random.Random(int(sys.argv[1])); N = int(sys.argv[2])
tot = bad = 0; hist = {}; nonface3 = 0; core3 = 0
for it in range(N):
    n = rng.randint(5, 7)
    g = planted(n, rng.randint(2, 5), rng, big=False)
    if g is None: continue
    A, H = g
    if soft_triangles(n, H): continue
    o = owner_map(H)
    inH = {v for t in H for v in t}
    anchors = [v for v in inH if not any(A[v, x] and frozenset((v, x)) not in o for x in range(n) if x != v)]
    if len(anchors) < 2: continue
    anchors = anchors[:4]
    sig, _ = sigma_bruteforce(A, H)
    # faces: x singleton in G xor sigma
    def singles(sg):
        B = A.copy()
        for (a, b, c), s in zip(H, sg):
            x, y = [(a, b), (b, c), (a, c)][s]; B[x, y] = B[y, x] = not B[x, y]
        return {v for v in range(n) if not any(B[v, x] for x in range(n) if x != v)}
    sing = [singles(s) for s in sig]
    face = lambda S: any(set(S) <= s for s in sing)
    for r in range(0, len(anchors) + 1):
        for S in itertools.combinations(anchors, r):
            if n + len(S) > 10: continue
            m = n + len(S)
            B = np.zeros((m, m), bool); B[:n, :n] = A
            for k, x in enumerate(S): B[x, n + k] = B[n + k, x] = True
            sg2, mn = sigma_bruteforce(B, H)
            actual = len(sg2) > 0
            pred = bool(sig) and all(face(c) for c in itertools.chain(itertools.combinations(S, 1), itertools.combinations(S, 2)))
            truth_face = face(S)
            tot += 1
            if actual != truth_face: bad += 1; print("LEMMA tinactive(a) mismatch", A.astype(int).tolist(), H, S)
            if actual != pred:
                bad += 1; print("A3 MISMATCH", A.astype(int).tolist(), H, S, actual, pred)
            if r >= 3 and not truth_face: nonface3 += 1
            if r >= 3 and not truth_face and pred: core3 += 1
print("checked (instance,S) pairs", tot, "mismatches", bad, "non-faces of size>=3", nonface3, "of which pairwise-consistent (3-core)", core3)

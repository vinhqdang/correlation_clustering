"""Independent classification of the binary relations R_AB (two members sharing exactly one vertex q, four hard cross pairs).
Built on real 5-vertex graphs; constraint = 'triple (p,q,r) is not an induced P3 after editing', evaluated on the edited GRAPH (no x/eps formula).
Also checks the full ZE set of the 5-vertex graph by partition enumeration, closure under all 729 rules, shapes."""
import itertools
import numpy as np
from collections import Counter
from a10_lib import *

def edited(A, H, s):
    B = A.copy()
    for (a, b, c), st in zip(H, s):
        x, y = [(a, b), (b, c), (a, c)][st]
        B[x, y] = B[y, x] = not B[x, y]
    return B

def p3(B, t):
    x, y, z = t
    return int(B[x, y]) + int(B[y, z]) + int(B[x, z]) == 2

rels = {}      # frozenset relation -> list of config descriptors
fullsig = {}
nconf = 0
for oA in itertools.permutations([0, 1, 2]):
    for oB in itertools.permutations([0, 3, 4]):
        for eps in itertools.product((0, 1), repeat=4):
            A = np.zeros((5, 5), bool)
            for (a, b, c) in (oA, oB):
                A[a, b] = A[b, a] = A[b, c] = A[c, b] = True
            for (p, r), e in zip([(1, 3), (1, 4), (2, 3), (2, 4)], eps):
                A[p, r] = A[r, p] = bool(e)
            H = [oA, oB]
            assert valid(A, H)
            nconf += 1
            R = set()
            for sA in range(3):
                for sB in range(3):
                    B = edited(A, H, (sA, sB))
                    if not any(p3(B, (p, 0, r)) for p in (1, 2) for r in (3, 4)):
                        R.add((sA, sB))
            rels.setdefault(frozenset(R), []).append((oA.index(0), oB.index(0), eps))
            sg, _ = sigma_bruteforce(A, H)
            # the 5-vertex Sigma = binary part  AND  unary parts (triples with one soft pair); check Sigma subset of R
            assert sg <= R
            fullsig[frozenset(sg)] = fullsig.get(frozenset(sg), 0) + 1
print("configs", nconf, "distinct binary relations", len(rels))
print("size histogram of relations:", sorted(Counter(len(r) for r in rels).items()))
print("distinct full 5-vertex Sigma:", len(fullsig), sorted(Counter(len(r) for r in fullsig).items()))
# all 729 majority rules explicitly
vals = [0, 1, 2]
perms = list(itertools.permutations(vals))
rules = []
for out in itertools.product(range(3), repeat=6):
    tbl = dict(zip(perms, out))
    def m(x, y, z, tbl=tbl):
        if x == y or x == z: return x
        if y == z: return y
        return tbl[(x, y, z)]
    rules.append(m)
print("rules", len(rules))
def closed_rel(R):
    R = list(R)
    for m in rules:
        for u, v, w in itertools.product(R, repeat=3):
            if (m(u[0], v[0], w[0]), m(u[1], v[1], w[1])) not in R: return False
    return True
bad = [R for R in rels if not closed_rel(R)]
print("relations not closed under some rule:", len(bad))
badf = [R for R in fullsig if not closed_rel(R)]
print("full sigma not closed:", len(badf))
# shape classification
def shape(R):
    n = len(R)
    rows = Counter(s for s, t in R); cols = Counter(t for s, t in R)
    if n <= 2: return f"size{n}"
    if n == 3:
        if max(rows.values()) == 3: return "line-row"
        if max(cols.values()) == 3: return "line-col"
        if all(s < 2 and t < 2 for s, t in R): return "3 in LR-square"
        return "other3:" + str(sorted(R))
    if n == 5:
        for p in range(3):
            for q in range(3):
                if R == {(p, t) for t in range(3)} | {(s, q) for s in range(3)}: return f"cross(p={p},q={q})"
        return "other5:" + str(sorted(R))
    return "other%d:" % n + str(sorted(R))
sh = Counter(shape(R) for R in rels)
for k, v in sorted(sh.items()): print(k, v)
# does a cross with p or q == 2 (M) appear?
for R in rels:
    if len(R)==3 and shape(R).startswith("line"): print(shape(R), sorted(R))

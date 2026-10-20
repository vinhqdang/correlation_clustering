"""a13_t8_b1.py: Theorem B1 with the csp-variant F' (member-adjacency locality), re-implemented from scratch, against exact OPT (subset DP).
Random graphs n<=11, random (not nec. maximal) modification-disjoint induced-P3 packings, budgets |H|+l, l=0..2, random deletion orders."""
import random, itertools, sys
from a13_lib import *
random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
def rand_inst(n):
    p = random.choice([.25, .35, .5])
    E = set(frozenset(q) for q in itertools.combinations(range(n), 2) if random.random() < p)
    A = adjm(n, E)
    P3 = [(a, b, c) for b in range(n) for a, c in itertools.combinations(range(n), 2) if a != b and c != b and A[a][b] and A[b][c] and not A[a][c]]
    random.shuffle(P3); used = set(); mem = []
    cap = random.choice([2, 3, 4, 99])
    for (a, b, c) in P3:
        ps = [frozenset((a, b)), frozenset((b, c)), frozenset((a, c))]
        if any(p in used for p in ps): continue
        used |= set(ps); mem.append((a, b, c))
        if len(mem) >= cap: break
    return dict(n=n, E=E, members=mem, T=[])
def solve(inst, l, depth=0):
    H = len(inst['members']); B = H + l
    n = inst['n']; cnt = counted_pairs(inst)
    far = [frozenset(p) for p in itertools.combinations(range(n), 2) if frozenset(p) not in cnt]
    f = Fprime(inst, far)
    if f is not None and f <= B: return True
    if l == 0: return False
    T = list(far); random.shuffle(T)
    for e in list(T):
        T2 = [x for x in T if x != e]
        g = Fprime(inst, T2)
        if g is None or g > B: T = T2
    for e in T:
        E2 = set(inst['E']) ^ {e}
        if solve(dict(inst, E=E2), l-1, depth+1): return True
    return False
tot = 0; mism = 0; seen = {}
for it in range(int(sys.argv[2]) if len(sys.argv) > 2 else 60):
    n = random.randint(6, 10); inst = rand_inst(n)
    H = len(inst['members'])
    if H == 0: continue
    opt = opt_exact(inst)
    for l in range(0, 3):
        r = solve(inst, l); truth = opt <= H + l
        tot += 1
        if r != truth: mism += 1; print('MISMATCH', it, l, opt, H)
        seen[(truth)] = seen.get(truth, 0) + 1
print('B1(F\') vs exact OPT: %d (instance,budget) pairs, %d mismatches; YES/NO split %s' % (tot, mism, seen))

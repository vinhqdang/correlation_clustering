"""w11_kappa.py -- unit-lock Helly number kappa_1 of every enumerated component relation (max size of a minimal inconsistent family of singleton locks)."""
import pickle, itertools
from collections import Counter
rels = pickle.load(open('/tmp/claude-0/-home-user-correlation-clustering/404aaf8c-18f1-5f19-bc11-1c2160450d1f/scratchpad/w11_rels.pkl', 'rb'))
def kappa(R):
    R = list(R)
    if not R: return 0
    k = len(R[0]); best = 0
    # partial assignments: for each coordinate value in {None,0,1,2}
    for pa in itertools.product((None, 0, 1, 2), repeat=k):
        locks = [(i, v) for i, v in enumerate(pa) if v is not None]
        if not locks: continue
        feas = any(all(t[i] == v for i, v in locks) for t in R)
        if feas: continue
        # minimal: every proper subset feasible
        minimal = True
        for j in range(len(locks)):
            sub = locks[:j] + locks[j + 1:]
            if not sub: continue
            if not any(all(t[i] == v for i, v in sub) for t in R): minimal = False; break
        if minimal: best = max(best, len(locks))
    return best
for name, rs in rels.items():
    c = Counter(kappa(R) for R in rs)
    print(name, len(rs), 'distribution of kappa_1:', dict(sorted(c.items())))

# sub-relations R ∩ prod(lists): lists any nonempty subset of {0,1,2}
subsets = [s for r in (1, 2, 3) for s in itertools.combinations((0, 1, 2), r)]
for name in ('k4', 'dia', 'k3', 'link'):
    worst = Counter(); seenR = set()
    for R in rels[name]:
        k = len(next(iter(R))) if R else 0
        for lists in itertools.product(subsets, repeat=k):
            Rp = frozenset(t for t in R if all(t[i] in lists[i] for i in range(k)))
            if Rp in seenR or not Rp: continue
            seenR.add(Rp)
            worst[kappa(Rp)] += 1
    print(name, 'sub-relations R∩lists: distinct', len(seenR), 'kappa_1 distribution', dict(sorted(worst.items())), flush=True)

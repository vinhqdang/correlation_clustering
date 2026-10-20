"""Is the set Sigma of zero-excess member-state vectors of a tight instance closed under a coordinatewise majority operation on {L,R,M}?
(Closure under a majority would imply 2-decomposability, hence lock cores of size <= 2.)  For each instance, all 729 majority rules are tested
(m(a,a,b)=a, and an arbitrary value for each of the 6 orderings of three distinct states).  Reports per c = max members per vertex."""
import sys, itertools, random, time, numpy as np, networkx as nx
from w6c_lib import *
from w6c_nflib import packings_of_size
PERMS = list(itertools.permutations(range(3)))
RULES = list(itertools.product(range(3), repeat=6))
def states(lab, H):
    out = []
    for (a, b, c) in H:
        t = (lab[a] == lab[b], lab[b] == lab[c], lab[a] == lab[c])
        out.append({(False, True, False): 0, (True, False, False): 1, (True, True, True): 2}[tuple(bool(x) for x in t)])
    return tuple(out)
def closed_rules(Sig):
    Sig = list(Sig); X = np.array(Sig)            # (s, p)
    s, p = X.shape
    if s < 3: return set(range(len(RULES)))
    good = set()
    # precompute for triples: coordinatewise; for speed use all ordered triples of distinct vectors only where some coordinate has all-distinct values
    idx = np.array(list(itertools.product(range(s), repeat=3)))
    A, B, C = X[idx[:, 0]], X[idx[:, 1]], X[idx[:, 2]]
    # coordinatewise majority where two agree
    maj = np.where(A == B, A, np.where(A == C, A, np.where(B == C, B, -1)))
    sigset = {tuple(v) for v in Sig}
    # rows with no -1: check closure independent of rule
    full = (maj >= 0).all(1)
    if full.any():
        for row in maj[full]:
            if tuple(row) not in sigset: return set()
    # rows with -1 entries depend on the rule
    rows = np.nonzero(~full)[0]
    if len(rows) == 0: return set(range(len(RULES)))
    # for each rule check
    pi = {pm: i for i, pm in enumerate(PERMS)}
    cand = []
    for ri, rule in enumerate(RULES):
        ok = True
        for r in rows:
            v = maj[r].copy()
            for q in np.nonzero(v < 0)[0]:
                v[q] = rule[pi[(A[r, q], B[r, q], C[r, q])]]
            if tuple(v) not in sigset: ok = False; break
        if ok: cand.append(ri)
    return set(cand)
if __name__ == '__main__':
    n = int(sys.argv[1]); T = float(sys.argv[2]); seed = int(sys.argv[3]); rng = random.Random(seed); t0 = time.time()
    stat = {}; ex = {}
    while time.time() - t0 < T:
        g = nx.gnm_random_graph(n, rng.randint(n, int(2.2 * n)), seed=rng.randrange(10 ** 9))
        if not nx.is_connected(g): continue
        A = np.array(nx.to_numpy_array(g), dtype=np.int8)
        opt, lab, same, cost, P = brute_opt(A)
        Hs = packings_of_size(A, opt)
        if not Hs: continue
        optidx = np.nonzero(cost == opt)[0]
        rng.shuffle(Hs)
        for H in Hs[:3]:
            cnt = [0] * n
            for t in H:
                for v in t: cnt[v] += 1
            c = max(cnt)
            Sig = {states(lab[i], H) for i in optidx}
            if len(Sig) > 40 or len(Sig) < 3: continue
            cr = closed_rules(Sig)
            st = stat.setdefault(c, [0, 0, set(range(len(RULES)))]); st[0] += 1
            if cr: st[1] += 1; st[2] &= cr
            else:
                if c not in ex or len(H) < len(ex[c][1]): ex[c] = (sorted(g.edges()), H, len(Sig))
    for c, (tot, ok, common) in sorted(stat.items()):
        print('n', n, 'c', c, 'instances', tot, 'closed under some majority', ok, 'rules good for all those:', len(common))
    for c, v in ex.items(): print('NOT majority-closed, c =', c, v)

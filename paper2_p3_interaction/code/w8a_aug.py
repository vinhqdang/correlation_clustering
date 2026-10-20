"""r-exchange augmentation local search on packings: remove <= r members, add >= r+1 pair-disjoint induced P3s on the freed pairs.  Compare local optimum with nu (ILP) and OPT."""
import random, sys, itertools, collections, time
from w8a_lib import *
from w8a_ilp import rand_graph
def best_fill(cands, pp, need):
    """exists >= need pair-disjoint P3s among cands (pp[t] = pair list)?  DFS"""
    cands = list(cands)
    def rec(i, used, cnt):
        if cnt >= need: return True
        if cnt + (len(cands) - i) < need: return False
        for j in range(i, len(cands)):
            t = cands[j]
            if not any(p in used for p in pp[t]):
                for p in pp[t]: used.add(p)
                ok = rec(j + 1, used, cnt + 1)
                for p in pp[t]: used.discard(p)
                if ok: return True
        return False
    return rec(0, set(), 0)
def fill(cands, pp, need):
    cands = list(cands); best = []
    def rec(i, used, cur):
        nonlocal best
        if len(cur) >= need and not best: best = list(cur); return
        if best: return
        for j in range(i, len(cands)):
            t = cands[j]
            if not any(p in used for p in pp[t]):
                for p in pp[t]: used.add(p)
                cur.append(t); rec(j + 1, used, cur); cur.pop()
                for p in pp[t]: used.discard(p)
                if best: return
    rec(0, set(), [])
    return best
def improve(g, H, r):
    P = g.p3s(); pp = {t: mp(t) for t in P}
    soft = soft_pairs(H)
    for k in range(1, r + 1):
        for R in itertools.combinations(range(len(H)), k):
            Hrest = [H[i] for i in range(len(H)) if i not in R]
            used = set(p for t in Hrest for p in pp.get(t, mp(t)))
            cands = [t for t in P if not any(p in used for p in pp[t])]
            add = fill(cands, pp, k + 1)
            if add: return Hrest + add
    return None
def local_opt(g, H, r):
    H = list(H); steps = 0
    while True:
        H2 = improve(g, H, r)
        if H2 is None: return H, steps
        H = maximalize(g, H2); steps += 1
if __name__ == '__main__':
    seed = int(sys.argv[1]); N = int(sys.argv[2]); rng = random.Random(seed); st = collections.Counter(); t0 = time.time()
    for it in range(N):
        n = rng.randint(9, 15); mode = rng.choice(['gnp', 'tree+chords', 'cycle+chords', 'cubic', 'planted'])
        g = G8(n, rand_graph(n, mode, rng))
        gg = nx.Graph(list(g.E))
        if len(g.E) < 3 or gg.number_of_nodes() < n or not nx.is_connected(gg): continue
        P = g.p3s()[:]; rng.shuffle(P); H = []; used = set()
        for t in P:
            if not any(p in used for p in mp(t)): H.append(t); used.update(mp(t))
        nu, _ = nu_mod(g); opt, _ = opt_ilp(g)
        res = []
        for r in (1, 2):
            Hr, steps = local_opt(g, H, r); res.append(len(Hr))
        st[('nu-|H|=%d' % (nu - len(H)), 'closed by r=1: %s' % (res[0] == nu), 'closed by r=2: %s' % (res[1] == nu))] += 1
    print('seed', seed, 'time %.0fs' % (time.time() - t0))
    for k, v in sorted(st.items(), key=str): print('  ', k, v)

"""Audit 11: a REAL implementation of the Theorem-PL dynamic program (promise-FPT DP of notes/wave8_repack.md Sec. 4), written from the
text of the proof (not from w8a_pl.py), plus brute-force ground truth over all set partitions.

DP: tree decomposition (T,beta) of I(H); bag state = (x_P in 7 nonempty subsets of pairs(P) for P in beta_t ;
E_v subset of Elig(v), |E_v|<=l, for every covered v in X_t cap L_B).  Checks at a bag: (S) symmetry, and every vertex triple whose three
pair statuses are visible at the bag.  Cost charged at the top bag of each item.  Value m* = min over root states.
Elig(v) = hard partners q in Ball(v) (radius 2 in G[L_B]) with (q covered, in L_B, co-bagged with v) or (q uncovered local hub).
Brute force: CLASS-OPT = min cost over ALL partitions whose hard edits are all eligible and |E_v|<=l for covered v.
"""
import sys, random, itertools
from a11_lib import *

class Inst:
    def __init__(s, g, H, l, bags, tree_edges, loose=False):
        s.g, s.H, s.l = g, H, l
        n = g.n
        s.own = soft_owner(H)
        s.cov = {v for t in H for v in t}
        cnt = [0] * n
        for t in H:
            for v in t: cnt[v] += 1
        s.c = max(cnt) if H else 0
        s.D = 2 * s.c
        s.B = 2 * s.D + 3 * l
        s.dh = [sum(1 for w in g.adj[v] if K(v, w) not in s.own) for v in range(n)]
        s.LB = {v for v in range(n) if s.dh[v] <= s.B}
        # Ball radius 2 in G[LB]
        s.ball = {}
        for v in s.LB:
            b = {v}; fr = {v}
            for _ in range(2):
                fr = {w for x in fr for w in g.adj[x] if w in s.LB} - b
                b |= fr
            s.ball[v] = b
        s.bags = bags            # list of sets of member indices
        s.tree = tree_edges      # list of (u,v)
        s.X = [set(v for P in bag for v in H[P]) for bag in bags]
        s.cob = set()
        for X in s.X:
            for a, b in itertools.combinations(sorted(X), 2): s.cob.add((a, b))
        def cobag(a, b): return K(a, b) in s.cob
        s.cobag = cobag
        s.hub = {}
        for r in range(n):
            if r in s.cov or r not in s.LB: continue
            vs = sorted(v for v in s.ball[r] if v in s.cov)
            s.hub[r] = all(cobag(a, b) for a, b in itertools.combinations(vs, 2))
        s.elig = {}
        for v in s.cov:
            if v not in s.LB: s.elig[v] = set(); continue
            el = set()
            for q in s.ball[v]:
                if q == v or K(v, q) in s.own: continue
                if q in s.cov:
                    if q in s.LB and (loose or cobag(v, q)): el.add(q)
                elif s.hub.get(q) or (loose and q in s.LB): el.add(q)
            s.elig[v] = el

    def hard_pair_class(s, p, q):
        """'cc' covered-covered eligible, 'hub' (r,v) eligible, None not editable"""
        if K(p, q) in s.own: return None
        if p in s.cov and q in s.cov: return 'cc' if q in s.elig[p] else None
        if p in s.cov and q not in s.cov: return 'hub' if q in s.elig[p] else None
        if q in s.cov and p not in s.cov: return 'hub' if p in s.elig[q] else None
        return None

SUB = [frozenset(c) for r in (1, 2, 3) for c in itertools.combinations(range(3), r)]   # 7 nonempty subsets of pair indices (0:ab,1:bc,2:ac)

def member_pairs(t):
    a, b, c = t
    return [K(a, b), K(b, c), K(a, c)]

def dp(I, check_S=True, restrict_label=True, verbose=False, coverage=True):
    g, H, l = I.g, I.H, I.l
    nb = len(I.bags)
    adjT = {i: [] for i in range(nb)}
    for u, v in I.tree: adjT[u].append(v); adjT[v].append(u)
    root = 0; parent = {root: None}; order = [root]
    for u in order:
        for v in adjT[u]:
            if v not in parent: parent[v] = u; order.append(v)
    children = {i: [v for v in adjT[i] if parent.get(v) == i] for i in range(nb)}
    n = g.n
    triples = list(itertools.combinations(range(n), 3))
    # coverage (visibility) check, state independent
    def lab_vertices(t): return sorted(v for v in I.X[t] if v in I.cov and v in I.LB)
    def visible(t, p, q):
        if K(p, q) in I.own: return I.own[K(p, q)] in I.bags[t]
        cl = I.hard_pair_class(p, q)
        if cl is None: return True
        if cl == 'cc': return p in I.X[t] or q in I.X[t]
        r, v = (p, q) if p not in I.cov else (q, p)
        return v in I.X[t]
    cov_fail = 0
    for T in triples:
        a, b, c = T
        if not any(all(visible(t, p, q) for p, q in ((a, b), (a, c), (b, c))) for t in range(nb)): cov_fail += 1
    tables = {}
    for t in reversed(order):
        mem = sorted(I.bags[t]); vs = lab_vertices(t)
        par = parent[t]
        new_mem = [P for P in mem if par is None or P not in I.bags[par]]
        new_vs = [v for v in vs if par is None or v not in I.X[par]]
        # labels per vertex: subsets of Elig of size <= l
        def label_opts(v):
            el = sorted(I.elig[v]); out = []
            for r in range(0, l + 1):
                for cmb in itertools.combinations(el, r): out.append(frozenset(cmb))
            return out
        lopts = [label_opts(v) for v in vs]
        # statuses helper
        def edited(p, q, x, E):
            """x: dict member->frozenset idx ; E: dict vertex->frozenset ; None if not visible"""
            key = K(p, q)
            if key in I.own:
                P = I.own[key]
                if P not in x: return None
                return member_pairs(H[P]).index(key) in x[P]
            cl = I.hard_pair_class(p, q)
            if cl is None: return False
            if cl == 'cc':
                if p in E: return q in E[p]
                if q in E: return p in E[q]
                return None
            r, v = (p, q) if p not in I.cov else (q, p)
            if v in E: return r in E[v]
            return None
        # children projections
        child_proj = []
        for ch in children[t]:
            cm = sorted(set(mem) & I.bags[ch]); cv = sorted(set(vs) & set(lab_vertices(ch)))
            proj = {}
            for (xs, es), cost in tables[ch].items():
                chmem = sorted(I.bags[ch]); chvs = lab_vertices(ch)
                key = (tuple(xs[chmem.index(P)] for P in cm), tuple(es[chvs.index(v)] for v in cv))
                if key not in proj or cost < proj[key]: proj[key] = cost
            child_proj.append((cm, cv, proj))
        table = {}
        for labs in itertools.product(*lopts):
            E = dict(zip(vs, labs))
            if check_S:
                bad = False
                for p, q in itertools.combinations(vs, 2):
                    if (q in E[p]) != (p in E[q]): bad = True; break
                if bad: continue
            lcost = sum(len(E[v] & I.cov) + 2 * len(E[v] - I.cov) for v in new_vs)
            # children label compatibility prefilter (labels part)
            for xs in itertools.product(SUB, repeat=len(mem)):
                x = dict(zip(mem, xs))
                ok = True
                for T in triples:
                    a, b, c = T
                    st = [edited(a, b, x, E), edited(a, c, x, E), edited(b, c, x, E)]
                    if None in st: continue
                    ne = sum(1 for (p, q), e in zip(((a, b), (a, c), (b, c)), st) if g.e(p, q) != e)
                    if ne == 2: ok = False; break
                if not ok: continue
                cost = lcost + sum(2 * len(x[P]) for P in new_mem)
                good = True
                for cm, cv, proj in child_proj:
                    key = (tuple(x[P] for P in cm), tuple(E[v] for v in cv))
                    if key not in proj: good = False; break
                    cost += proj[key]
                if good:
                    table[(tuple(x[P] for P in mem), tuple(E[v] for v in vs))] = cost
        tables[t] = table
    best = min(tables[root].values()) if tables[root] else None
    return (best // 2 if best is not None else None), cov_fail

def brute_class_opt(I):
    g, H, l = I.g, I.H, I.l
    D = diff_matrix(g)
    cost = D.sum(1)
    ok = np.ones(len(D), bool)
    pi = g.pidx
    cnt = np.zeros((len(D), g.n), dtype=np.int16)
    for (p, q), j in pi.items():
        if K(p, q) in I.own: continue
        cl = I.hard_pair_class(p, q)
        if cl is None:
            ok &= ~D[:, j]
        else:
            for v in (p, q):
                if v in I.cov: cnt[:, v] += D[:, j]
    ok &= (cnt <= l).all(1)
    if not ok.any(): return None, None
    return int(cost[ok].min()), ok, cost, D

def elimination_decomposition(Iadj, order):
    """tree decomposition from an elimination order of the graph Iadj (dict)"""
    adj = {v: set(Iadj[v]) for v in Iadj}
    pos = {v: i for i, v in enumerate(order)}
    bags = {}; par = {}
    for v in order:
        nbrs = set(adj[v])
        bags[v] = nbrs | {v}
        for a in nbrs:
            adj[a] |= nbrs - {a}
            adj[a].discard(v)
        for a in nbrs: pass
        par[v] = min(nbrs, key=lambda u: pos[u]) if nbrs else None
    idx = {v: i for i, v in enumerate(order)}
    bl = [bags[v] for v in order]
    edges = []
    roots = [v for v in order if par[v] is None]
    for v in order:
        if par[v] is not None: edges.append((idx[v], idx[par[v]]))
    for a, b in zip(roots, roots[1:]): edges.append((idx[a], idx[b]))   # join forest components
    return bl, edges

def run_random(seed, trials, maxw=2, verbose=True):
    rng = random.Random(seed)
    stats = dict(inst=0, sound_viol=0, dp_vs_class_mismatch=0, promise_true=0, promise_complete_fail=0, promise_false=0, cov_fail=0, strictly_worse=0, skipped=0)
    for _ in range(trials):
        n = rng.randint(6, 9); p = rng.choice([0.25, 0.35, 0.5])
        g = Gr(n, [(a, b) for a in range(n) for b in range(a + 1, n) if rng.random() < p])
        if not g.p3(): continue
        H = random_packing(g, rng)
        o, opts = opt_and_optima(g)
        if len(H) > 5: stats['skipped'] += 1; continue
        l = max(1, o - len(H)) if rng.random() < 0.7 else o - len(H)
        Ia = I_graph(H)
        order = list(Ia); rng.shuffle(order)
        bags, edges = elimination_decomposition(Ia, order)
        if max((len(b) for b in bags), default=0) - 1 > maxw: stats['skipped'] += 1; continue
        # decomposition may contain isolated empty case
        inst = Inst(g, H, l, bags, edges)
        val, cf = dp(inst)
        cb = brute_class_opt(inst)
        stats['inst'] += 1; stats['cov_fail'] += cf > 0
        cls = cb[0]
        if val is not None and val < o: stats['sound_viol'] += 1; print('SOUNDNESS VIOLATION', seed, sorted(g.E), H, l)
        if val != cls: stats['dp_vs_class_mismatch'] += 1; print('MISMATCH dp', val, 'class', cls, 'opt', o, sorted(g.E), H, l, bags)
        # promise: some optimum (size o <= |H|+l) in the class
        _, ok, cost, D = cb if cls is not None else (None, None, None, None)
        prom = cls is not None and cls == o
        if prom:
            stats['promise_true'] += 1
            if val != o: stats['promise_complete_fail'] += 1
        else:
            stats['promise_false'] += 1
            if val is None or val > o: stats['strictly_worse'] += 1
    print(seed, stats)



def sparse_instance(rng, n):
    """random tree + few chords"""
    ed = [(rng.randrange(i), i) for i in range(1, n)]
    for _ in range(rng.randint(1, 3)):
        a, b = rng.sample(range(n), 2); ed.append((a, b))
    return Gr(n, ed)

def run_sparse(seed, trials, loose=False, check_S=True):
    rng = random.Random(seed)
    st = dict(inst=0, sound_viol=0, mismatch=0, promise_true=0, promise_complete_fail=0, promise_false=0, cov_fail=0, nonclass_gap=0, skipped=0)
    for _ in range(trials):
        n = rng.randint(7, 9)
        g = sparse_instance(rng, n)
        if not g.p3(): continue
        H = random_packing(g, rng)
        if len(H) < 3 or len(H) > 6: st['skipped'] += 1; continue
        o, opts = opt_and_optima(g)
        l = max(1, o - len(H))
        if o - len(H) > 2: st['skipped'] += 1; continue
        Ia = I_graph(H)
        order = sorted(Ia, key=lambda P: min(H[P])) if rng.random() < 0.5 else rng.sample(list(Ia), len(Ia))
        bags, edges = elimination_decomposition(Ia, order)
        if max((len(b) for b in bags), default=0) - 1 > 2: st['skipped'] += 1; continue
        inst = Inst(g, H, l, bags, edges, loose=loose)
        val, cf = dp(inst, check_S=check_S)
        cb = brute_class_opt(inst)
        cls = cb[0]
        st['inst'] += 1; st['cov_fail'] += cf > 0
        if val is not None and val < o: st['sound_viol'] += 1
        if val != cls: st['mismatch'] += 1
        prom = cls == o
        if prom:
            st['promise_true'] += 1
            if val != o: st['promise_complete_fail'] += 1
        else:
            st['promise_false'] += 1
            if val is None or val > o: st['nonclass_gap'] += 1
    print('sparse', seed, 'loose' if loose else '', 'noS' if not check_S else '', st)

if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'rand': run_random(int(sys.argv[2]), int(sys.argv[3]))
    elif mode == 'sparse': run_sparse(int(sys.argv[2]), int(sys.argv[3]))
    elif mode == 'loose': run_sparse(int(sys.argv[2]), int(sys.argv[3]), loose=True)
    elif mode == 'noS': run_sparse(int(sys.argv[2]), int(sys.argv[3]), check_S=False)

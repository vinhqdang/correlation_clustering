"""Try to realise a state-level L/R 3-core by leaf lock members.  For a tight instance with a pairwise-consistent non-extendable L/R assignment (coords i0,i1,i2 ; values v)
attach to each core member Pk a leaf member Qk = (x_k, ., .) sharing one vertex with Pk, x_k a clean anchor.  For each option compute by ILP the set of states of Pk reachable with
x_k locked; keep options whose forced set is exactly {v_k}.  Then test the three anchors: pairwise faces but jointly not (lock 3-core)."""
import sys, itertools, random, time, numpy as np, networkx as nx
from w6c_lib import *
from w6c_nflib import packings_of_size
from w6c_majority import states
ST = {0: (0, 1, 0), 1: (1, 0, 0), 2: (1, 1, 1)}
def sfix(mem, s):
    a, b, c = mem; ab, bc, ac = ST[s]; return {key(a, b): ab, key(b, c): bc, key(a, c): ac}
def augment(A, H, members_opts):
    """members_opts: list of (k, v, variant, pattern); returns new adjacency, new packing, anchors"""
    n = len(A); E = {key(i, j) for i in range(n) for j in range(i + 1, n) if A[i][j]}
    Hn = list(H); nn = n; anchors = []
    for (k, v, var, pat) in members_opts:
        Pk = H[k]; u, x = nn, nn + 1; nn += 2
        others = [w for w in Pk if w != v]
        Q = [(x, u, v), (x, v, u), (v, u, x), (u, v, x)][var]
        Hn.append(Q); E |= {key(Q[0], Q[1]), key(Q[1], Q[2])}
        for idx, (w) in enumerate(others):
            for jj in range(1):
                pass
        bits = [(pat >> i) & 1 for i in range(2)]
        for w, b in zip(others, bits):
            if b: E.add(key(u, w))
        anchors.append(x)
    return adj_from_edges(nn, E), Hn, anchors
def allowed_states(A2, Hn, Pk, x):
    n2 = len(A2); out = set()
    lock = {key(x, t): 0 for t in range(n2) if t != x}
    for s in (0, 1, 2):
        fx = dict(lock); fx.update(sfix(Pk, s))
        r = ilp_opt(A2, fixed=fx)
        if r is not None and r[0] == len(Hn): out.add(s)
    return out
if __name__ == '__main__':
    n = int(sys.argv[1]); T = float(sys.argv[2]); seed = int(sys.argv[3]); rng = random.Random(seed); t0 = time.time()
    inst = 0; opts_ok = 0; found = 0
    while time.time() - t0 < T:
        g = nx.gnm_random_graph(n, rng.randint(n, int(2.3 * n)), seed=rng.randrange(10 ** 9))
        if not nx.is_connected(g): continue
        A = np.array(nx.to_numpy_array(g), dtype=np.int8)
        opt, lab, same, cost, P = brute_opt(A)
        Hs = packings_of_size(A, opt)
        if not Hs: continue
        optidx = np.nonzero(cost == opt)[0]
        H = rng.choice(Hs); p = len(H)
        if p < 3: continue
        Sig = {states(lab[i], H) for i in optidx}
        S = np.array(sorted(Sig)); core = None
        for coords in itertools.combinations(range(p), 3):
            sub = S[:, coords]; trip = {tuple(r) for r in sub}
            pairs = {(i, j): {(r[i], r[j]) for r in sub} for i, j in ((0, 1), (0, 2), (1, 2))}
            for v in itertools.product((0, 1), repeat=3):
                if v not in trip and all((v[i], v[j]) in pairs[(i, j)] for i, j in ((0, 1), (0, 2), (1, 2))): core = (coords, v); break
            if core: break
        if not core: continue
        inst += 1
        coords, vals = core
        good = []
        for k, vk in zip(coords, vals):
            gk = []
            for v in H[k]:
                for var in range(4):
                    for pat in range(4):
                        A2, Hn, anc = augment(A, H, [(k, v, var, pat)])
                        ok, _ = check_packing(A2, Hn)
                        if not ok or ilp_opt(A2)[0] != len(Hn): continue
                        al = allowed_states(A2, Hn, H[k], anc[0])
                        if al == {vk}: gk.append((k, v, var, pat))
            good.append(gk)
        print('instance', inst, 'core', core, 'options forcing the core state:', [len(x) for x in good], flush=True)
        if all(good):
            opts_ok += 1
            tried = 0
            for combo in itertools.product(*[gk[:6] for gk in good]):
                A3, H3, anc = augment(A, H, list(combo))
                ok, _ = check_packing(A3, H3)
                if not ok: continue
                o3 = ilp_opt(A3)[0]
                if o3 != len(H3): continue
                def fc(Ss):
                    fx = {}
                    for a in Ss:
                        for t in range(len(A3)):
                            if t != a: fx[key(a, t)] = 0
                    r = ilp_opt(A3, fixed=fx); return r is not None and r[0] == o3
                pf = [fc(list(c_)) for c_ in itertools.combinations(anc, 2)]
                if all(pf) and not fc(anc):
                    found += 1; print('LOCK 3-CORE!', sorted(g.edges()), H, combo, flush=True)
                tried += 1
            print('   combos tried', tried, flush=True)
    print('n', n, 'core instances', inst, 'with forcing options for all three members', opts_ok, 'lock 3-cores', found, flush=True)

"""w11_maj2.py -- constraint scopes of Sigma_0 for c<=2 as triangle-connected components of I(H): K3, diamond, K4, and links (edges in no triangle).
Closure of each distinct relation under the M-majority m(x,y,z)=x if two equal else M(=2), and under all 729 majorities."""
import itertools, sys, time
import numpy as np
from collections import Counter
TOG = np.array([[0,1,0],[1,0,0],[1,1,1]])   # TOG[state][k], k: 0=ab,1=bc,2=ac

def centre_variants(vs):
    return [(vs[(i + 1) % 3], vs[i], vs[(i + 2) % 3]) for i in range(3)]

def build(members_vs, nv, extra_X_pairs_allowed=()):
    """members_vs: list of vertex lists (3 each). Enumerate roles (centres) and hard bits; return dict relation->first config."""
    m = len(members_vs)
    soft = {}
    rels = {}
    allpairs = list(itertools.combinations(range(nv), 2))
    softpairs_vs = set()
    for vs in members_vs:
        for p in itertools.combinations(sorted(vs), 2): softpairs_vs.add(p)
    hardpairs = [p for p in allpairs if p not in softpairs_vs]
    triples = list(itertools.combinations(range(nv), 3))
    states = list(itertools.product((0,1,2), repeat=m))
    S = np.array(states)                           # (ns, m)
    for roles in itertools.product(*[centre_variants(vs) for vs in members_vs]):
        soft = {}
        for i, (a, b, c) in enumerate(roles):
            for p, k in (((a,b),0),((b,c),1),((a,c),2)):
                soft[tuple(sorted(p))] = (i, k)
        # per triple: soft contribution array and list of hard pairs
        tinfo = []
        for t in triples:
            ps = [tuple(sorted(q)) for q in itertools.combinations(t, 2)]
            contrib = np.zeros(len(states), dtype=np.int8)
            hp = []
            for p in ps:
                if p in soft:
                    i, k = soft[p]; contrib += TOG[S[:, i], k]
                else:
                    hp.append(p)
            tinfo.append((contrib, hp))
        # group triples by hard-pair tuple
        nh = len(hardpairs)
        hidx = {p: j for j, p in enumerate(hardpairs)}
        xopts = [(0, 1, 2)] if True else None
        # allowed arrays per triple per constant sum
        for bits in itertools.product((0, 1, 2) if extra_X_pairs_allowed else (0, 1), repeat=nh):
            # bit value 2 = X (only allowed for pairs listed)
            ok = True
            for j, p in enumerate(hardpairs):
                if bits[j] == 2 and p not in extra_X_pairs_allowed: ok = False; break
            if not ok: continue
            allowed = np.ones(len(states), dtype=bool)
            for contrib, hp in tinfo:
                if any(bits[hidx[p]] == 2 for p in hp): continue
                cs = contrib + sum(bits[hidx[p]] for p in hp)
                allowed &= (cs != 2)
            R = frozenset(map(tuple, S[allowed].tolist()))
            if R not in rels: rels[R] = (roles, bits)
    return rels

def maj_table(vals):
    perms = list(itertools.permutations((0, 1, 2)))
    tab = np.zeros((3, 3, 3), dtype=np.int8)
    for x, y, z in itertools.product((0, 1, 2), repeat=3):
        if x == y or x == z: tab[x, y, z] = x
        elif y == z: tab[x, y, z] = y
    for p, v in zip(perms, vals): tab[p] = v
    return tab

def closed(R, tab):
    A = np.array(sorted(R), dtype=np.int8)
    if len(A) == 0: return True
    n, m = A.shape
    code = lambda X: (X.astype(np.int64) * (3 ** np.arange(m))).sum(-1)
    inR = np.zeros(3 ** m, dtype=bool); inR[code(A)] = True
    idx = np.arange(n)
    I1, I2, I3 = np.meshgrid(idx, idx, idx, indexing='ij')
    out = np.stack([tab[A[I1.ravel(), c], A[I2.ravel(), c], A[I3.ravel(), c]] for c in range(m)], axis=-1)
    return bool(inR[code(out)].all())

def report(name, rels, Mtab, ops=None):
    t0 = time.time()
    bad = [R for R in rels if not closed(R, Mtab)]
    print(f'{name}: distinct relations {len(rels)}, not closed under M-majority: {len(bad)}  ({time.time()-t0:.0f}s)', flush=True)
    if ops is not None:
        good = [v for v, tab in ops if all(closed(R, tab) for R in rels)]
        print(f'   majorities (of 729) preserving all {name} relations: {len(good)}', flush=True)
    return bad

if __name__ == '__main__':
    Mvals = (2,) * 6
    Mtab = maj_table(Mvals)
    ops = [(v, maj_table(v)) for v in itertools.product((0, 1, 2), repeat=6)]
    # link: two members sharing one vertex, cross pairs hard bits only (no X)
    links = build([[0, 1, 2], [0, 3, 4]], 5)
    report('link', links, Mtab, ops)
    k3 = build([[0, 1, 3], [0, 2, 4], [1, 2, 5]], 6)
    report('K3', k3, Mtab, ops)
    # diamond: A=(0,1,2) B=(0,3,4) C=(1,3,5) D=(2,4,6); pair (5,6) may be X (owned by an outside member)
    t0 = time.time()
    dia = build([[0, 1, 2], [0, 3, 4], [1, 3, 5], [2, 4, 6]], 7, extra_X_pairs_allowed=((5, 6),))
    print('diamond relations enumerated', len(dia), round(time.time() - t0), 's', flush=True)
    bad = report('diamond', dia, Mtab, ops)
    k4 = build([[0, 1, 2], [0, 3, 4], [1, 3, 5], [2, 4, 5 + 1]], 7) if False else None
    # K4: vertices x_ij; members A=(AB,AC,AD) B=(AB,BC,BD) C=(AC,BC,CD) D=(AD,BD,CD): ids AB0 AC1 AD2 BC3 BD4 CD5
    k4 = build([[0, 1, 2], [0, 3, 4], [1, 3, 5], [2, 4, 5]], 6)
    report('K4', k4, Mtab, ops)
    import pickle
    pickle.dump({'link': links, 'k3': k3, 'dia': dia, 'k4': k4}, open('/tmp/claude-0/-home-user-correlation-clustering/404aaf8c-18f1-5f19-bc11-1c2160450d1f/scratchpad/w11_rels.pkl', 'wb'))

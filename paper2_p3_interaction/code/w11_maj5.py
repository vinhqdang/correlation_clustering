"""w11_maj5.py -- five-state (member-only) relations over D5={L,R,M,4,5} for the c<=2 components with hard pairs constant;
search (SAT/CEGAR) for an operation m on the 'sunflower patterns' (all-distinct triples with <=1 defective value, output non-defective)
such that for every relation R and every triple u,v,w in R^3 whose coordinates are each 'all equal or at most one defective',
m(u,v,w) in R.  (Weighted Helly gluing step.)"""
import itertools, sys, time, random
import numpy as np
from pysat.solvers import Cadical153
from w11_maj2 import centre_variants

# states: 0=L(ab),1=R(bc),2=M(ac),3=4 (ab,bc),4=5 (ab,bc,ac); togetherness (ab,bc,ac) after edit
TOG5 = np.array([[0,1,0],[1,0,0],[1,1,1],[0,0,0],[0,0,1]])
DEF = {3, 4}

def rel5(members, nv, hard, X=()):
    m = len(members)
    soft = {}
    for i, (a, b, c) in enumerate(members):
        for p, k in (((a,b),0),((b,c),1),((a,c),2)): soft[tuple(sorted(p))] = (i, k)
    states = np.array(list(itertools.product(range(5), repeat=m)))
    allowed = np.ones(len(states), dtype=bool)
    for t in itertools.combinations(range(nv), 3):
        ps = [tuple(sorted(q)) for q in itertools.combinations(t, 2)]
        if any(p in X for p in ps): continue
        cs = np.zeros(len(states), dtype=np.int8)
        for p in ps:
            if p in soft:
                i, k = soft[p]; cs += TOG5[states[:, i], k]
            else: cs += hard[p]
        allowed &= (cs != 2)
    return states[allowed]

def configs_link(maxn=None):
    out = []
    cross = [(1,3),(1,4),(2,3),(2,4)]
    for ma in centre_variants([0,1,2]):
        for mb in centre_variants([0,3,4]):
            for bits in itertools.product((0,1), repeat=4):
                hard = {tuple(sorted(p)): b for p, b in zip(cross, bits)}
                out.append(([ma, mb], 5, hard, ()))
    return out

def configs_k3():
    out = []
    A=[0,1,3]; B=[0,2,4]; C=[1,2,5]
    soft = set()
    for vs in (A,B,C):
        for p in itertools.combinations(sorted(vs),2): soft.add(p)
    hp = [p for p in itertools.combinations(range(6),2) if p not in soft]
    for roles in itertools.product(*[centre_variants(vs) for vs in (A,B,C)]):
        for bits in itertools.product((0,1), repeat=6):
            out.append((list(roles), 6, dict(zip(hp,bits)), ()))
    return out

def configs_diamond(rnd, count):
    out = []
    mem = [[0,1,2],[0,3,4],[1,3,5],[2,4,6]]
    soft = set()
    for vs in mem:
        for p in itertools.combinations(sorted(vs),2): soft.add(p)
    hp = [p for p in itertools.combinations(range(7),2) if p not in soft]
    for _ in range(count):
        roles = [rnd.choice(centre_variants(vs)) for vs in mem]
        bits = [rnd.choice((0,1)) for _ in hp]
        X = ()
        hard = dict(zip(hp, bits))
        if rnd.random() < 0.3: X = ((5,6),)
        out.append((roles, 7, hard, X))
    return out

def configs_k4(rnd, count):
    out = []
    mem = [[0,1,2],[0,3,4],[1,3,5],[2,4,5]]
    soft = set()
    for vs in mem:
        for p in itertools.combinations(sorted(vs),2): soft.add(p)
    hp = [p for p in itertools.combinations(range(6),2) if p not in soft]
    for roles in itertools.product(*[centre_variants(vs) for vs in mem]):
        for bits in itertools.product((0,1), repeat=len(hp)):
            out.append((list(roles), 6, dict(zip(hp,bits)), ()))
    return out

NO = int(sys.argv[1]) if len(sys.argv) > 1 else 3
PAT_ID = {}
def pat_vars():
    vals = [0,1,2,3,4]
    pats = []
    for t in itertools.permutations(vals, 3):
        nd = sum(1 for x in t if x in DEF)
        if nd <= 1: pats.append(t)
    return pats

def check_relation(R, mtab, limit=3_000_000, rnd=None):
    """return list of violating (pattern tuple lists) for relation R (array nxk) ; mtab dict pattern->value. Triples restricted to sunflower-type."""
    n, k = R.shape
    viol = []
    idx = np.arange(n)
    # full cube if small, else random sample
    if n ** 3 <= limit:
        I1, I2, I3 = np.meshgrid(idx, idx, idx, indexing='ij'); I1 = I1.ravel(); I2 = I2.ravel(); I3 = I3.ravel()
    else:
        I1 = np.random.randint(0, n, limit); I2 = np.random.randint(0, n, limit); I3 = np.random.randint(0, n, limit)
    U = R[I1]; V = R[I2]; W = R[I3]          # (T,k)
    ok = np.ones(len(I1), dtype=bool)
    defmask = lambda X: (X == 3) | (X == 4)
    for c in range(k):
        eq = (U[:, c] == V[:, c]) & (V[:, c] == W[:, c])
        nd = defmask(U[:, c]).astype(int) + defmask(V[:, c]) + defmask(W[:, c])
        # allowed: all equal, or <=1 defective  (if two equal-but-not-all with defect count 2 -> excluded)
        ok &= eq | (nd <= 1)
    return U[ok], V[ok], W[ok]

if __name__ == '__main__':
    rnd = random.Random(1); np.random.seed(1)
    pats = pat_vars()
    print('patterns', len(pats))
    solver = Cadical153()
    var = {}
    for p in pats:
        for o in range(NO):
            var[(p, o)] = len(var) + 1
        solver.add_clause([var[(p, o)] for o in range(NO)])
        for o1 in range(NO):
            for o2 in range(o1 + 1, NO): solver.add_clause([-var[(p, o1)], -var[(p, o2)]])
    def base_out(a, b, c):
        if a == b or a == c: return a
        if b == c: return b
        return None
    rels = []
    t0 = time.time()
    for name, cf in (('link', configs_link()), ('k3', configs_k3())):
        seen = set()
        for (mem, nv, hard, X) in cf:
            R = rel5(mem, nv, hard, X)
            key = R.tobytes() + bytes([R.shape[1]])
            if key in seen or len(R) == 0: continue
            seen.add(key); rels.append((name, R))
        print(name, 'distinct D5 relations', len(seen), round(time.time() - t0), 's', flush=True)
    dia = configs_diamond(rnd, 400)
    seen = set()
    for (mem, nv, hard, X) in dia:
        R = rel5(mem, nv, hard, X)
        key = R.tobytes()
        if key in seen or len(R) == 0: continue
        seen.add(key); rels.append(('dia', R))
    print('diamond sampled distinct', len(seen), flush=True)
    sizes = [len(R) for _, R in rels]
    print('relations', len(rels), 'max size', max(sizes))
    # CEGAR
    cur = {}
    it = 0
    while True:
        it += 1
        if not solver.solve():
            print('UNSAT: no operation on sunflower patterns works (after', it, 'iterations)'); break
        model = set(l for l in solver.get_model() if l > 0)
        mt = {}
        for p in pats:
            for o in range(NO):
                if var[(p, o)] in model: mt[p] = o
        def mfun(a, b, c):
            r = base_out(a, b, c)
            return r if r is not None else mt[(a, b, c)]
        newc = 0
        for name, R in rels:
            U, V, W = check_relation(R, mt)
            if len(U) == 0: continue
            k = R.shape[1]
            Rset = set(map(tuple, R.tolist()))
            # evaluate m on all triples
            bad = 0
            for u, v, w in zip(U.tolist(), V.tolist(), W.tolist()):
                out = tuple(mfun(u[c], v[c], w[c]) for c in range(k))
                if out not in Rset:
                    # clause: forbid the current values of the pattern vars involved
                    lits = []
                    for c in range(k):
                        p = (u[c], v[c], w[c])
                        if base_out(*p) is None: lits.append(-var[(p, mt[p])])
                    solver.add_clause(lits)
                    newc += 1; bad += 1
                    if bad > 200: break
        print('iter', it, 'new clauses', newc, round(time.time() - t0), 's', flush=True)
        if newc == 0:
            print('FOUND operation m (on sunflower patterns) closing all tested D5 relations'); print(mt); break

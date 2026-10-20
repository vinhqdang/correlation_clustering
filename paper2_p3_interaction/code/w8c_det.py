"""Wave 8: SAT/CEGIS synthesis of a *detector* gadget on a plain P3 spine.

Spine: members P_j=(s_{2j}, s_{2j+1}, s_{2j+2}), j=0..q-1 (states L,R,M), an induced path.  Valid ZE state vectors of the bare
spine are the factors R^a L^b and R^a M L^b.  Defect group at position j (1<=j<=q-1):
    F_j = {sigma: sigma_{j-1}=R and sigma_j=L}  U  {sigma: sigma_j=M}.
Goal: extra vertices W (anchor x = first extra vertex, clean, end of exactly one member) and extra members such that
  (U) every valid sigma extends to a ZE solution of the whole instance (anchor unconstrained);
  (L) every valid sigma NOT in F_j extends with x a singleton; every valid sigma in F_j does NOT.
usage: python3 w8c_det.py q j w Mx cmax [seconds] [seed]   (Mx = number of extra members)
"""
import sys, os, time, random, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w8c_lib import *
from pysat.solvers import Cadical153
from pysat.card import CardEnc, EncType


def valid_sigmas(q):
    out = []
    for a in range(q + 1):
        out.append(tuple(['R'] * a + ['L'] * (q - a)))
    for a in range(q):
        out.append(tuple(['R'] * a + ['M'] + ['L'] * (q - a - 1)))
    return out


FMODE = int(os.environ.get('FMODE', '0'))


def in_F(sig, j):
    if FMODE == 0:
        return (sig[j - 1] == 'R' and sig[j] == 'L') or sig[j] == 'M'
    if FMODE == 1:   # sanity: lock forbids sigma_j = L
        return sig[j] == 'L'
    if FMODE == 2:   # sanity: lock forbids sigma_j = M
        return sig[j] == 'M'
    if FMODE == 3:   # forbid only the flip at j (allow M at j)
        return sig[j - 1] == 'R' and sig[j] == 'L'
    if FMODE == 4:   # forbid only M at j
        return sig[j] == 'M'
    if FMODE == 5:
        return False
    if FMODE == 6:
        return True
    if FMODE == 7:   # lock forces flip at j
        return not (sig[j - 1] == 'R' and sig[j] == 'L')



STATE = {'L': (0, 1, 0), 'R': (1, 0, 0), 'M': (1, 1, 1)}


class Synth:
    def __init__(self, q, j, w, Mx, cmax, extra_edges_to_spine=True):
        self.q, self.j, self.w, self.Mx, self.cmax = q, j, w, Mx, cmax
        ns = 2 * q + 1
        self.ns = ns
        self.n = n = ns + w
        self.x = ns  # anchor
        self.spine_members = [(2 * i, 2 * i + 1, 2 * i + 2) for i in range(q)]
        vp = self.vp = VarPool()
        self.TRUE = vp.new()
        cl = self.cl = [[self.TRUE]]
        P = self.P = list(itertools.combinations(range(n), 2))
        spine_soft = set()
        for m in self.spine_members:
            spine_soft |= pairs_of(m)
        spine_edges = set(key(i, i + 1) for i in range(ns - 1))
        self.E = {}
        for p in P:
            if p[1] < ns:  # spine-spine pair: fixed
                self.E[p] = self.TRUE if p in spine_edges else -self.TRUE
            else:
                self.E[p] = vp.new()
        # extra-member candidates
        self.Mem = {}
        for b in range(n):
            if b == self.x:
                continue
            for a, c in itertools.combinations([v for v in range(n) if v != b], 2):
                # all three pairs must not be spine-soft pairs; a pair with both ends in spine may only be ac (E fixed 0)
                ps = (key(a, b), key(b, c), key(a, c))
                if any(p in spine_soft for p in ps):
                    continue
                if key(a, b)[1] < ns and key(a, b) not in spine_edges:
                    continue  # ab must be an edge
                if key(b, c)[1] < ns and key(b, c) not in spine_edges:
                    continue
                if key(a, c)[1] < ns and key(a, c) in spine_edges:
                    continue
                # a pair with both ends spine and edge: not allowed unless ... it's a soft spine pair: excluded above
                self.Mem[(b, a, c)] = vp.new()
        cover = {p: [] for p in P}
        inmem = {v: [] for v in range(n)}
        for (b, a, c), m in self.Mem.items():
            cl += [[-m, self.E[key(a, b)]], [-m, self.E[key(b, c)]], [-m, -self.E[key(a, c)]]]
            for p in (key(a, b), key(b, c), key(a, c)):
                cover[p].append(m)
            for v in (a, b, c):
                inmem[v].append(m)
        self.S = {}
        for p in P:
            if p in spine_soft:
                self.S[p] = self.TRUE
                continue
            s = vp.new()
            self.S[p] = s
            cl.append([-s] + cover[p])
            for m in cover[p]:
                cl.append([-m, s])
            for m1, m2 in itertools.combinations(cover[p], 2):
                cl.append([-m1, -m2])
        self.H0, self.H1 = {}, {}
        for p in P:
            h0, h1 = vp.new(), vp.new()
            self.H0[p], self.H1[p] = h0, h1
            S_, E_ = self.S[p], self.E[p]
            cl += [[-h0, -S_], [-h0, -E_], [h0, S_, E_]]
            cl += [[-h1, -S_], [-h1, E_], [h1, S_, -E_]]
        allm = list(self.Mem.values())
        enc = CardEnc.equals(lits=allm, bound=Mx, top_id=vp.n, encoding=EncType.seqcounter)
        vp.n = max(vp.n, enc.nv)
        cl += enc.clauses
        # c bound: count spine members too
        spine_cnt = {v: sum(1 for m in self.spine_members if v in m) for v in range(n)}
        for v in range(n):
            bnd = cmax - spine_cnt[v]
            if bnd < 0:
                bnd = 0
            enc = CardEnc.atmost(lits=inmem[v], bound=bnd, top_id=vp.n, encoding=EncType.seqcounter)
            vp.n = max(vp.n, enc.nv)
            cl += enc.clauses
        # anchor: in exactly one member, as end; clean
        x = self.x
        enc = CardEnc.equals(lits=inmem[x], bound=1, top_id=vp.n, encoding=EncType.seqcounter)
        vp.n = max(vp.n, enc.nv)
        cl += enc.clauses
        for v in range(n):
            if v != x:
                cl.append([-self.E[key(x, v)], self.S[key(x, v)]])
        # scenarios: existential
        self.sigmas = valid_sigmas(q)
        for sig in self.sigmas:
            self.add_scenario(sig, singleton=False)
            if not in_F(sig, j):
                self.add_scenario(sig, singleton=True)

    def add_scenario(self, sig, singleton):
        n, vp, cl = self.n, self.vp, self.cl
        T = {}
        for p in self.P:
            T[p] = vp.new()
        for u, v, w_ in itertools.combinations(range(n), 3):
            a, b, c = T[(u, v)], T[(v, w_)], T[(u, w_)]
            cl += [[-a, -b, c], [-a, -c, b], [-b, -c, a]]
        for p in T:
            cl += [[self.S[p], -T[p], self.E[p]], [self.S[p], T[p], -self.E[p]]]
        for (b, a, c), m in self.Mem.items():
            e1, e2, e3 = -T[key(a, b)], -T[key(b, c)], T[key(a, c)]
            cl += [[-m, e1, e2, e3], [-m, -e1, -e2], [-m, -e1, -e3], [-m, -e2, -e3]]
        for i, st in enumerate(sig):
            a, b, c = self.spine_members[i]
            tab, tbc, tac = STATE[st]
            for (p, val) in ((key(a, b), tab), (key(b, c), tbc), (key(a, c), tac)):
                cl.append([T[p]] if val else [-T[p]])
        if singleton:
            for v in range(n):
                if v != self.x:
                    cl.append([-T[key(self.x, v)]])

    def decode(self, model):
        ms = set(l for l in model if l > 0)
        E = set()
        for p, v in self.E.items():
            if p[1] < self.ns:
                if v == self.TRUE:
                    E.add(p)
            elif v in ms:
                E.add(p)
        members = list(self.spine_members)
        for (b, a, c), v in self.Mem.items():
            if v in ms:
                members.append((a, b, c))
        return E, members

    def verify(self, E, members):
        """return (sig, lab) counterexample: sigma in F with locked extension, or None"""
        for sig in self.sigmas:
            if not in_F(sig, self.j):
                continue
            fa = []
            for i, st in enumerate(sig):
                a, b, c = self.spine_members[i]
                tab, tbc, tac = STATE[st]
                for (p, val) in ((key(a, b), tab), (key(b, c), tbc), (key(a, c), tac)):
                    (fa if not val else fa)
            fixed_t = []
            fixed_a = []
            for i, st in enumerate(sig):
                a, b, c = self.spine_members[i]
                tab, tbc, tac = STATE[st]
                for (p, val) in ((key(a, b), tab), (key(b, c), tbc), (key(a, c), tac)):
                    (fixed_t if val else fixed_a).append(p)
            lab = ze_solve(self.n, E, members, singletons=[self.x], fixed_apart=fixed_a, fixed_together=fixed_t)
            if lab is not None:
                return sig, lab
        return None

    def cex_clause(self, lab, perm=None):
        n = self.n
        lab2 = lab
        if perm is not None:
            lab2 = [None] * n
            for v in range(n):
                lab2[perm[v]] = lab[v]
        same = lambda u, v: lab2[u] == lab2[v]
        c = []
        for p in self.P:
            if p[1] < self.ns:
                continue
            c.append(self.H0[p] if same(*p) else self.H1[p])
        for (b, a, cc), m in self.Mem.items():
            cnt = (0 if same(a, b) else 1) + (0 if same(b, cc) else 1) + (1 if same(a, cc) else 0)
            if cnt != 1:
                c.append(m)
        # spine-spine pairs that could become soft/hard: H0/H1 constant w.r.t. spine fixed; spine pairs that are non-soft are
        # consistent by construction of sigma, so omitted.
        return c


def run(q, j, w, Mx, cmax, seconds, seed):
    syn = Synth(q, j, w, Mx, cmax)
    rng = random.Random(seed)
    solver = Cadical153(bootstrap_with=syn.cl)
    t0 = time.time()
    it = 0
    extras = list(range(syn.ns + 1, syn.n))
    while time.time() - t0 < seconds:
        if not solver.solve():
            print('UNSAT after', it, 'cex: no detector in this model (q=%d j=%d w=%d Mx=%d cmax=%d)' % (q, j, w, Mx, cmax), flush=True)
            return None
        E, members = syn.decode(solver.get_model())
        r = syn.verify(E, members)
        it += 1
        if r is None:
            print('FOUND after', it, flush=True)
            return E, members, syn
        sig, lab = r
        for rep in range(30):
            perm = list(range(syn.n))
            if rep:
                e2 = extras[:]
                rng.shuffle(e2)
                for a_, b_ in zip(extras, e2):
                    perm[a_] = b_
            # spine relabelling is not a symmetry (spine fixed); only extra vertices are permuted
            solver.add_clause(syn.cex_clause(lab, perm))
        if it % 50 == 0:
            print('iter', it, round(time.time() - t0, 1), flush=True)
    print('TIMEOUT', it, flush=True)
    return 'timeout'


if __name__ == '__main__':
    q, j, w, Mx, cmax = map(int, sys.argv[1:6])
    sec = float(sys.argv[6]) if len(sys.argv) > 6 else 600
    seed = int(sys.argv[7]) if len(sys.argv) > 7 else 1
    res = run(q, j, w, Mx, cmax, sec, seed)
    if res and res != 'timeout':
        E, members, syn = res
        print('n', syn.n, 'E', sorted(E))
        print('members', members)
        print('valid', check_instance(syn.n, E, members))
        for sig in syn.sigmas:
            fa, ft = [], []
            for i, st in enumerate(sig):
                a, b, c = syn.spine_members[i]
                tab, tbc, tac = STATE[st]
                for (p, val) in ((key(a, b), tab), (key(b, c), tbc), (key(a, c), tac)):
                    (ft if val else fa).append(p)
            u = ze_solve(syn.n, E, members, fixed_apart=fa, fixed_together=ft) is not None
            l = ze_solve(syn.n, E, members, singletons=[syn.x], fixed_apart=fa, fixed_together=ft) is not None
            print(''.join(sig), 'unlocked', u, 'locked', l, 'inF', in_F(sig, syn.j))

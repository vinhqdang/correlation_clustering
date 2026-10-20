"""Wave 8: SAT/CEGIS synthesis of a *comparator* between two one-hot chains OH_m (A and B).
Spec (threshold equality at position t):  ZE solutions of A+B+gadget have spine-state pair (a,b) iff [a<=t] == [b<=t].
Vertices: A (cells), B, W extra.  New members may use any vertices (pairs inside one chain that are base-soft are excluded).
usage: python3 w8c_cmp.py m t w Mx cmax [seconds] [seed]      env PURE=1 (no hard edges at new pairs)  WIN=1 (cells t-WIN..t+WIN+1)
"""
import sys, os, time, random, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w8c_lib import *
from w8c_oh import oh_instance, spine_state_vector
from pysat.solvers import Cadical153
from pysat.card import CardEnc, EncType

STATE = {'L': (0, 1, 0), 'R': (1, 0, 0), 'M': (1, 1, 1)}


def cell_of(name):
    return name[1] // 2 if name[0] == 's' else name[1]


class CmpSynth:
    def __init__(self, m, t, w, Mx, cmax, win=1, pure=False, mode='cmp'):
        n0, E0, mem0, names0, idx0 = oh_instance(m)
        self.m, self.t, self.mode = m, t, mode
        self.n0 = n0
        # vertex ids: chain A: 0..n0-1, chain B: n0..2n0-1, extras 2n0..
        self.nA = n0
        n = self.n = 2 * n0 + w
        self.w = w
        self.base_members = [tuple(v for v in tr) for tr in mem0] + [tuple(v + n0 for v in tr) for tr in mem0]
        self.base_E = set(E0) | set((u + n0, v + n0) for (u, v) in E0)
        self.names0, self.idx0 = names0, idx0
        self.base_soft = set()
        for tr in self.base_members:
            self.base_soft |= pairs_of(tr)
        # active window cells
        lo, hi = t - win, t + win + 1
        self.active = set()
        for v in range(n0):
            if lo <= cell_of(names0[v]) <= hi:
                self.active.add(v)
                self.active.add(v + n0)
        for x in range(2 * n0, n):
            self.active.add(x)
        vp = self.vp = VarPool()
        self.TRUE = vp.new()
        cl = self.cl = [[self.TRUE]]
        P = self.P = list(itertools.combinations(range(n), 2))
        base_vs = set(range(2 * n0))
        self.E = {}
        self.fixed = {}
        for p in P:
            u, v = p
            same_chain = (u < n0 and v < n0) or (u >= n0 and v >= n0 and v < 2 * n0)
            if u < 2 * n0 and v < 2 * n0 and same_chain:
                self.E[p] = self.TRUE if p in self.base_E else -self.TRUE
                self.fixed[p] = True
            elif not (u in self.active and v in self.active):
                self.E[p] = -self.TRUE
                self.fixed[p] = True
            else:
                self.E[p] = vp.new()
        self.Mem = {}
        for b in range(n):
            if b not in self.active:
                continue
            for a_, c_ in itertools.combinations([v for v in self.active if v != b], 2):
                ps = (key(a_, b), key(b, c_), key(a_, c_))
                if any(q in self.base_soft for q in ps):
                    continue
                if any(self.fixed.get(q) for q in ps[:2]) and any((q in self.fixed and self.E[q] == -self.TRUE) for q in ps[:2]):
                    continue
                if (ps[2] in self.fixed and self.E[ps[2]] == self.TRUE):
                    continue
                if (ps[0] in self.fixed and self.E[ps[0]] == -self.TRUE) or (ps[1] in self.fixed and self.E[ps[1]] == -self.TRUE):
                    continue
                # at least one pair must be non-base (new): guaranteed since base-soft excluded and same-chain base pairs are fixed
                self.Mem[(b, a_, c_)] = vp.new()
        cover = {p: [] for p in P}
        inmem = {v: [] for v in range(n)}
        for (b, a_, c_), mv in self.Mem.items():
            cl += [[-mv, self.E[key(a_, b)]], [-mv, self.E[key(b, c_)]], [-mv, -self.E[key(a_, c_)]]]
            for q in (key(a_, b), key(b, c_), key(a_, c_)):
                cover[q].append(mv)
            for v in (a_, b, c_):
                inmem[v].append(mv)
        self.S = {}
        for p in P:
            if p in self.base_soft:
                self.S[p] = self.TRUE
                continue
            s = vp.new()
            self.S[p] = s
            cl.append([-s] + cover[p])
            for mv in cover[p]:
                cl.append([-mv, s])
            for m1, m2 in itertools.combinations(cover[p], 2):
                cl.append([-m1, -m2])
        if pure:
            for p in P:
                if p not in self.fixed and p not in self.base_soft:
                    cl.append([self.S[p], -self.E[p]])
        self.H0, self.H1 = {}, {}
        for p in P:
            h0, h1 = vp.new(), vp.new()
            self.H0[p], self.H1[p] = h0, h1
            S_, E_ = self.S[p], self.E[p]
            cl += [[-h0, -S_], [-h0, -E_], [h0, S_, E_]]
            cl += [[-h1, -S_], [-h1, E_], [h1, S_, -E_]]
        allm = list(self.Mem.values())
        enc = CardEnc.atmost(lits=allm, bound=Mx, top_id=vp.n, encoding=EncType.seqcounter)
        vp.n = max(vp.n, enc.nv)
        cl += enc.clauses
        cnt0 = {v: 0 for v in range(n)}
        for tr in self.base_members:
            for v in tr:
                cnt0[v] += 1
        for v in range(n):
            bnd = max(0, cmax - cnt0[v])
            enc = CardEnc.atmost(lits=list(dict.fromkeys(inmem[v])), bound=bnd, top_id=vp.n, encoding=EncType.seqcounter)
            vp.n = max(vp.n, enc.nv)
            cl += enc.clauses
        self.spineA = [(idx0[('s', 2 * j)], idx0[('s', 2 * j + 1)], idx0[('s', 2 * j + 2)]) for j in range(m)]
        self.spineB = [tuple(v + n0 for v in tr) for tr in self.spineA]
        # scenarios
        self.must, self.mustnot = [], []
        for a in range(m):
            for b in range(m):
                ok = ((a <= t) == (b <= t)) if mode == 'cmp' else True
                (self.must if ok else self.mustnot).append((a, b))
        for (a, b) in self.must:
            self.add_scenario(a, b)

    def fixlists(self, a, b):
        m = self.m
        sa = ''.join(['R'] * a + ['M'] + ['L'] * (m - a - 1))
        sb = ''.join(['R'] * b + ['M'] + ['L'] * (m - b - 1))
        fa, ft = [], []
        for spine, s in ((self.spineA, sa), (self.spineB, sb)):
            for tr, st in zip(spine, s):
                x, y, z = tr
                tab, tbc, tac = STATE[st]
                for (p, val) in ((key(x, y), tab), (key(y, z), tbc), (key(x, z), tac)):
                    (ft if val else fa).append(p)
        return fa, ft

    def add_scenario(self, a, b):
        n, vp, cl = self.n, self.vp, self.cl
        T = {p: vp.new() for p in self.P}
        for u, v, w_ in itertools.combinations(range(n), 3):
            x, y, z = T[(u, v)], T[(v, w_)], T[(u, w_)]
            cl += [[-x, -y, z], [-x, -z, y], [-y, -z, x]]
        for p in T:
            cl += [[self.S[p], -T[p], self.E[p]], [self.S[p], T[p], -self.E[p]]]
        for (b_, a_, c_), mv in self.Mem.items():
            e1, e2, e3 = -T[key(a_, b_)], -T[key(b_, c_)], T[key(a_, c_)]
            cl += [[-mv, e1, e2, e3], [-mv, -e1, -e2], [-mv, -e1, -e3], [-mv, -e2, -e3]]
        for tr in self.base_members:
            x, y, z = tr
            e1, e2, e3 = -T[key(x, y)], -T[key(y, z)], T[key(x, z)]
            cl += [[e1, e2, e3], [-e1, -e2], [-e1, -e3], [-e2, -e3]]
        fa, ft = self.fixlists(a, b)
        for p in ft:
            cl.append([T[p]])
        for p in fa:
            cl.append([-T[p]])

    def decode(self, model):
        ms = set(l for l in model if l > 0)
        E = set()
        for p in self.P:
            e = self.E[p]
            if e == self.TRUE or (e > 0 and e in ms):
                E.add(p)
        members = list(self.base_members)
        for (b, a_, c_), mv in self.Mem.items():
            if mv in ms:
                members.append((a_, b, c_))
        return E, members

    def cex(self, E, members):
        for (a, b) in self.mustnot:
            fa, ft = self.fixlists(a, b)
            lab = ze_solve(self.n, E, members, fixed_apart=fa, fixed_together=ft)
            if lab is not None:
                c = []
                same = lambda u, v: lab[u] == lab[v]
                for p in self.P:
                    if p in self.fixed:
                        continue
                    c.append(self.H0[p] if same(*p) else self.H1[p])
                for (b_, a_, c_), mv in self.Mem.items():
                    cnt = (0 if same(a_, b_) else 1) + (0 if same(b_, c_) else 1) + (1 if same(a_, c_) else 0)
                    if cnt != 1:
                        c.append(mv)
                return c
        return None


def run(m, t, w, Mx, cmax, seconds, seed, win, pure, mode='cmp'):
    syn = CmpSynth(m, t, w, Mx, cmax, win=win, pure=pure, mode=mode)
    print('built: n', syn.n, 'candidate members', len(syn.Mem), 'vars', syn.vp.n, 'clauses', len(syn.cl), flush=True)
    solver = Cadical153(bootstrap_with=syn.cl)
    t0 = time.time()
    it = 0
    while time.time() - t0 < seconds:
        if not solver.solve():
            print('UNSAT after', it, 'cex', flush=True)
            return None
        E, members = syn.decode(solver.get_model())
        it += 1
        c = syn.cex(E, members)
        if c is None:
            print('FOUND after', it, flush=True)
            return syn, E, members
        solver.add_clause(c)
        if it % 50 == 0:
            print('iter', it, round(time.time() - t0, 1), flush=True)
    print('TIMEOUT', it, flush=True)
    return 'timeout'


if __name__ == '__main__':
    m, t, w, Mx, cmax = map(int, sys.argv[1:6])
    sec = float(sys.argv[6]) if len(sys.argv) > 6 else 600
    seed = int(sys.argv[7]) if len(sys.argv) > 7 else 1
    res = run(m, t, w, Mx, cmax, sec, seed, int(os.environ.get('WIN', '1')), os.environ.get('PURE', '0') == '1', os.environ.get('MODE', 'cmp'))
    if res and res != 'timeout':
        syn, E, members = res
        print('valid', check_instance(syn.n, E, members))
        newm = members[len(syn.base_members):]
        print('new members', newm)
        print('new hard/soft edges among active:', sorted(p for p in E if p not in syn.fixed or False)[:80])
        # verify full table
        table = {}
        for a in range(syn.m):
            for b in range(syn.m):
                fa, ft = syn.fixlists(a, b)
                table[(a, b)] = ze_solve(syn.n, E, members, fixed_apart=fa, fixed_together=ft) is not None
        print('relation table (rows a, cols b):')
        for a in range(syn.m):
            print(''.join('1' if table[(a, b)] else '.' for b in range(syn.m)))

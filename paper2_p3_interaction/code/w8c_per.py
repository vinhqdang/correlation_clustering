"""Wave 8: periodic SAT/CEGIS synthesis on a skeleton chain.

Skeleton: m spine members P_j=(s_{2j}, s_{2j+1}, s_{2j+2}), vertex (cell j, class 0) = s_{2j}, (cell j, class 1) = s_{2j+1},
cell m has only class 0 (and twin).  Further classes: 2 = hard twin of class 0, 3 = hard twin of class 1 (optional glue),
classes >=4 are free extra vertices, class XC (if anchors) is the clean anchor of cell j (j<m).
Unknown: hard-pair pattern (translation invariant in the middle region, free near the two ends) and the extra members.
Specifications (SPEC):
  'pfg'   : all anchors locked -> infeasible;  all but cell i locked -> feasible for every interior i in RANGE
  'mspine': thick spine only: spine state vectors that extend == valid sigma without adjacent RL
usage: python3 w8c_per.py SPEC m cls cmax Mmem seconds [seed]
   cls = number of classes (>=2), anchors: env ANCH=1 puts an anchor class as the last class.
"""
import sys, os, time, random, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w8c_lib import *
from pysat.solvers import Cadical153
from pysat.card import CardEnc, EncType

STATE = {'L': (0, 1, 0), 'R': (1, 0, 0), 'M': (1, 1, 1)}


def valid_sigmas(q):
    out = []
    for a in range(q + 1):
        out.append(tuple(['R'] * a + ['L'] * (q - a)))
    for a in range(q):
        out.append(tuple(['R'] * a + ['M'] + ['L'] * (q - a - 1)))
    return out


def has_RL(sig):
    return any(sig[i] == 'R' and sig[i + 1] == 'L' for i in range(len(sig) - 1))


class Shared:
    def __init__(self):
        self.vp = VarPool()
        self.TRUE = self.vp.new()
        self.cl = [[self.TRUE]]
        self.Etype = {}
        self.Mkey = {}


class PerSynth:
    def __init__(self, m, ncls, cmax, Mmem, anchors=False, glue=True, lreg=2, rreg=2, span=2, extra_member_limit=None, shared=None, base=None):
        self.m, self.ncls, self.cmax, self.Mmem = m, ncls, cmax, Mmem
        self.base = base
        self.nospine = os.environ.get('NOSPINE', '0') == '1'
        self.anchors = anchors
        self.XC = ncls - 1 if anchors else None
        self.glue = glue
        self.lreg, self.rreg, self.span = lreg, rreg, span
        if shared is None:
            shared = Shared()
        self.shared = shared
        vp = self.vp = shared.vp
        self.TRUE = shared.TRUE
        cl = self.cl = shared.cl
        # vertices
        self.vid = {}
        self.vname = []
        for j in range(m + 1):
            for c in range(ncls):
                if self.exists(j, c):
                    self.vid[(j, c)] = len(self.vname)
                    self.vname.append((j, c))
        self.n = n = len(self.vname)
        P = self.P = list(itertools.combinations(range(n), 2))
        # spine members
        self.spine_members = []
        self.nospine = os.environ.get('NOSPINE', '0') == '1'
        if base is None and not self.nospine:
            for j in range(m):
                self.spine_members.append((self.vid[(j, 0)], self.vid[(j, 1)], self.vid[(j + 1, 0)]))
        elif base is not None:
            for (a_, b_, c_) in base['members']:
                self.spine_members.append((self.vid[a_], self.vid[b_], self.vid[c_]))
        self.spine_soft = set()
        for mm in self.spine_members:
            self.spine_soft |= pairs_of(mm)
        # fixed pair values
        self.fixedE = {}
        self.base_pairs = set()
        if base is not None:
            nb = base['nbase']
            for p in P:
                u_, v_ = self.vname[p[0]], self.vname[p[1]]
                if u_[1] < nb and v_[1] < nb:
                    self.base_pairs.add(p)
                    self.fixedE[p] = frozenset((u_, v_)) in base['E']
        for p in P:
            if self.nospine and base is None:
                (j, a), (k, b) = self.vname[p[0]], self.vname[p[1]]
                if k - j > span:
                    self.fixedE[p] = False
                continue
            if base is not None:
                if p in self.base_pairs:
                    continue
                (j, a), (k, b) = self.vname[p[0]], self.vname[p[1]]
                if k - j > span:
                    self.fixedE[p] = False
                continue
            (j, a), (k, b) = self.vname[p[0]], self.vname[p[1]]
            if k - j > span:
                self.fixedE[p] = False
                continue
            if p in self.spine_soft:
                self.fixedE[p] = True if (key(*p) in set(key(m_[0], m_[1]) for m_ in self.spine_members) | set(key(m_[1], m_[2]) for m_ in self.spine_members)) else False
                continue
            if a in (0, 1) and b in (0, 1):
                self.fixedE[p] = False  # spine-spine non-member pair: hard non-edge
                continue
            if glue and ((a, b) == (0, 2) or (a, b) == (2, 0) or (a, b) == (1, 3) or (a, b) == (3, 1)) and j == k:
                self.fixedE[p] = True
                self.glue_pairs = getattr(self, 'glue_pairs', set()) | {p}
                continue
        self.glue_pairs = getattr(self, 'glue_pairs', set())
        # edge variables
        self.Etype = shared.Etype
        self.E = {}
        for p in P:
            if p in self.fixedE:
                self.E[p] = self.TRUE if self.fixedE[p] else -self.TRUE
            else:
                tk = self.pkey(p)
                if tk not in self.Etype:
                    self.Etype[tk] = vp.new()
                self.E[p] = self.Etype[tk]
        # member candidates
        self.Mem = {}      # instance -> var
        self.Mkey = shared.Mkey
        for b in range(n):
            if anchors and self.vname[b][1] == self.XC:
                continue
            for a, c in itertools.combinations([v for v in range(n) if v != b], 2):
                ps = (key(a, b), key(b, c), key(a, c))
                if any(q in self.spine_soft for q in ps):
                    continue
                if any(q in self.glue_pairs for q in ps):
                    continue
                if any(q in self.base_pairs for q in ps):
                    continue
                if self.fixedE.get(ps[0]) is False or self.fixedE.get(ps[1]) is False:
                    continue
                if ps[2] in self.fixedE and self.fixedE[ps[2]] is True:
                    continue
                cells = [self.vname[v][0] for v in (a, b, c)]
                if max(cells) - min(cells) > span:
                    continue
                mk = self.mkey((b, a, c))
                if mk not in self.Mkey:
                    self.Mkey[mk] = vp.new()
                self.Mem[(b, a, c)] = self.Mkey[mk]
        cover = {p: [] for p in P}
        inmem = {v: [] for v in range(n)}
        centre_of = {v: [] for v in range(n)}
        for (b, a, c), mv in self.Mem.items():
            cl += [[-mv, self.E[key(a, b)]], [-mv, self.E[key(b, c)]], [-mv, -self.E[key(a, c)]]]
            for p in (key(a, b), key(b, c), key(a, c)):
                cover[p].append(mv)
            for v in (a, b, c):
                inmem[v].append(mv)
            centre_of[b].append(mv)
        self.S = {}
        for p in P:
            if p in self.spine_soft:
                self.S[p] = self.TRUE
                continue
            cv = list(dict.fromkeys(cover[p]))
            s = vp.new()
            self.S[p] = s
            cl.append([-s] + cv)
            for mv in cv:
                cl.append([-mv, s])
            # at most one covering member instance (distinct instances); shared var may appear twice for same instance? no.
            for m1, m2 in itertools.combinations(cv, 2):
                cl.append([-m1, -m2])
        if os.environ.get('PURE', '0') == '1':
            for p in P:
                if p not in self.fixedE and p not in self.spine_soft:
                    cl.append([self.S[p], -self.E[p]])      # no hard edges: every edge is soft
        self.H0, self.H1 = {}, {}
        for p in P:
            h0, h1 = vp.new(), vp.new()
            self.H0[p], self.H1[p] = h0, h1
            S_, E_ = self.S[p], self.E[p]
            cl += [[-h0, -S_], [-h0, -E_], [h0, S_, E_]]
            cl += [[-h1, -S_], [-h1, E_], [h1, S_, -E_]]
        # cardinality of members: total number of distinct member *instances* <= Mmem * (cells) is awkward; bound per-vertex c and total
        allm = list(dict.fromkeys(self.Mem.values()))
        # bound number of true template variables in the middle by Mmem (rough control of size)
        enc = CardEnc.atmost(lits=allm, bound=Mmem, top_id=vp.n, encoding=EncType.seqcounter)
        vp.n = max(vp.n, enc.nv)
        cl += enc.clauses
        spine_cnt = {v: 0 for v in range(n)}
        for mm in self.spine_members:
            for v in mm:
                spine_cnt[v] += 1
        for v in range(n):
            bnd = max(0, cmax - spine_cnt[v])
            lits = list(dict.fromkeys(inmem[v]))
            enc = CardEnc.atmost(lits=lits, bound=bnd, top_id=vp.n, encoding=EncType.seqcounter)
            vp.n = max(vp.n, enc.nv)
            cl += enc.clauses
        # anchors: exactly one member (as end), clean
        self.anchor_v = []
        if anchors:
            for j in range(m):
                x = self.vid[(j, self.XC)]
                self.anchor_v.append(x)
                lits = list(dict.fromkeys(inmem[x]))
                enc = CardEnc.equals(lits=lits, bound=1, top_id=vp.n, encoding=EncType.seqcounter)
                vp.n = max(vp.n, enc.nv)
                cl += enc.clauses
                for v in range(n):
                    if v != x:
                        cl.append([-self.E[key(x, v)], self.S[key(x, v)]])

    def exists(self, j, c):
        if self.base is not None:
            if c < self.base['nbase']:
                return not (j == self.m and c in (1, 3))
            return j < self.m
        if self.nospine:
            return not (self.anchors and c == self.XC and j == self.m)
        if j == self.m and c in (1, 3):
            return False
        if self.anchors and c == self.XC and j == self.m:
            return False
        if c == 3 and not self.glue and False:
            return False
        return True

    def region(self, j):
        if j < self.lreg:
            return ('L', j)
        if j > self.m - self.rreg:
            return ('R', j)
        return ('M',)

    def pkey(self, p):
        (j, a), (k, b) = self.vname[p[0]], self.vname[p[1]]
        rj, rk = self.region(j), self.region(k)
        if rj == ('M',) and rk == ('M',):
            return ('M', a, b, k - j)
        if rj[0] == 'R' or rk[0] == 'R':
            return ('R', j - self.m, a, k - self.m, b)
        return ('L', j, a, k, b)

    def mkey(self, mem):
        b, a, c = mem
        vs = [self.vname[v] for v in (b, a, c)]
        cells = [x[0] for x in vs]
        if all(self.region(j) == ('M',) for j in cells):
            lo = min(cells)
            return ('M', tuple((x[0] - lo, x[1]) for x in vs))
        if any(self.region(j)[0] == 'R' for j in cells):
            return ('R', tuple((x[0] - self.m, x[1]) for x in vs))
        return ('L', tuple(vs))

    # ---- scenarios
    def add_scenario(self, sig=None, locked=(), spine_free=False):
        vp, cl, n = self.vp, self.cl, self.n
        T = {p: vp.new() for p in self.P}
        for u, v, w_ in itertools.combinations(range(n), 3):
            a, b, c = T[(u, v)], T[(v, w_)], T[(u, w_)]
            cl += [[-a, -b, c], [-a, -c, b], [-b, -c, a]]
        for p in T:
            cl += [[self.S[p], -T[p], self.E[p]], [self.S[p], T[p], -self.E[p]]]
        for (b, a, c), mv in self.Mem.items():
            e1, e2, e3 = -T[key(a, b)], -T[key(b, c)], T[key(a, c)]
            cl += [[-mv, e1, e2, e3], [-mv, -e1, -e2], [-mv, -e1, -e3], [-mv, -e2, -e3]]
        for (a, b, c) in self.spine_members:
            e1, e2, e3 = -T[key(a, b)], -T[key(b, c)], T[key(a, c)]
            cl += [[e1, e2, e3], [-e1, -e2], [-e1, -e3], [-e2, -e3]]
        if sig is not None:
            for i, st in enumerate(sig):
                a, b, c = self.spine_members[i]
                tab, tbc, tac = STATE[st]
                for (p, val) in ((key(a, b), tab), (key(b, c), tbc), (key(a, c), tac)):
                    cl.append([T[p]] if val else [-T[p]])
        for j in locked:
            x = self.anchor_v[j]
            for v in range(n):
                if v != x:
                    cl.append([-T[key(x, v)]])
        return T

    def decode(self, model):
        ms = set(l for l in model if l > 0)
        E = set()
        for p in self.P:
            e = self.E[p]
            if (e == self.TRUE) or (e > 0 and e in ms):
                E.add(p)
        members = list(self.spine_members)
        for (b, a, c), mv in self.Mem.items():
            if mv in ms:
                members.append((a, b, c))
        return E, members

    def cex_clause(self, lab):
        c = []
        same = lambda u, v: lab[u] == lab[v]
        for p in self.P:
            if p in self.fixedE:
                if (self.fixedE[p] and not same(*p)) or ((not self.fixedE[p]) and same(*p)):
                    if p not in self.spine_soft:
                        # hard fixed pair violated by lab: this partition is never consistent -> cannot happen for lab from solver
                        pass
                continue
            c.append(self.H0[p] if same(*p) else self.H1[p])
        for (b, a, cc), mv in self.Mem.items():
            cnt = (0 if same(a, b) else 1) + (0 if same(b, cc) else 1) + (1 if same(a, cc) else 0)
            if cnt != 1:
                c.append(mv)
        return c


def spine_states2(syn, lab):
    out = []
    for (a, b, c) in syn.spine_members:
        tab, tbc, tac = lab[a] == lab[b], lab[b] == lab[c], lab[a] == lab[c]
        out.append({(1, 0, 0): 'R', (0, 1, 0): 'L', (1, 1, 1): 'M'}.get((int(tab), int(tbc), int(tac)), '?'))
    return ''.join(out)


def run(spec, m, ncls, cmax, Mmem, seconds, seed, anchors, glue=True):
    syn = PerSynth(m, ncls, cmax, Mmem, anchors=anchors, glue=glue)
    rng = random.Random(seed)
    sigmas = valid_sigmas(m)
    if spec == 'mspine':
        for s in sigmas:
            if not has_RL(s):
                syn.add_scenario(sig=s)
    elif spec == 'mexact':
        for s in sigmas:
            if s.count('M') == 1:
                syn.add_scenario(sig=s)
    elif spec == 'pfg':
        for i in range(0, m):
            syn.add_scenario(locked=[k for k in range(m) if k != i])
    solver = Cadical153(bootstrap_with=syn.cl)
    t0 = time.time()
    it = 0
    while time.time() - t0 < seconds:
        if not solver.solve():
            print('UNSAT after', it, 'cex', flush=True)
            return None
        E, members = syn.decode(solver.get_model())
        it += 1
        found = None
        if spec in ('mspine', 'mexact'):
            for s in sigmas:
                if (has_RL(s) if spec == 'mspine' else s.count('M') != 1):
                    fa, ft = [], []
                    for i, st in enumerate(s):
                        a, b, c = syn.spine_members[i]
                        tab, tbc, tac = STATE[st]
                        for (p, val) in ((key(a, b), tab), (key(b, c), tbc), (key(a, c), tac)):
                            (ft if val else fa).append(p)
                    lab = ze_solve(syn.n, E, members, fixed_apart=fa, fixed_together=ft)
                    if lab is not None:
                        found = lab
                        break
        elif spec == 'pfg':
            lab = ze_solve(syn.n, E, members, singletons=syn.anchor_v)
            if lab is not None:
                found = lab
        if found is None:
            print('FOUND after', it, flush=True)
            syn.model_set = set(l for l in solver.get_model() if l > 0)
            return syn, E, members
        solver.add_clause(syn.cex_clause(found))
        if it % 100 == 0:
            print('iter', it, round(time.time() - t0, 1), flush=True)
    print('TIMEOUT', it, flush=True)
    return 'timeout'


def translate_key(tk, m_new, m_orig):
    """map a template key of the m_new-chain to the key of the m_orig-chain (see notes: left/right anchored)."""
    if tk[0] == 'M':
        return tk
    if tk[0] == 'B':
        # pair key ('B', j, a, k, b) or member key ('B', ((cell,class),...)) -- members use tuple
        if isinstance(tk[1], tuple):
            cells = [c for (c, _) in tk[1]]
            shift = 0 if min(cells) <= 1 else (m_new - m_orig)
            return ('B', tuple((c - shift, cl) for (c, cl) in tk[1]))
        _, j, a, k, b = tk
        shift = 0 if min(j, k) <= 1 else (m_new - m_orig)
        return ('B', j - shift, a, k - shift, b)
    return tk


def extend_check(syn, model_set, m_new, ncls, cmax, Mmem, glue, anchors):
    syn2 = PerSynth(m_new, ncls, cmax, Mmem, anchors=anchors, glue=glue)
    ms = model_set
    E = set()
    for p in syn2.P:
        if p in syn2.fixedE:
            if syn2.fixedE[p]:
                E.add(p)
            continue
        tk = translate_key(syn2.pkey(p), m_new, syn.m)
        v = syn.Etype.get(tk)
        if v is None:
            raise KeyError(('pair', tk))
        if v in ms:
            E.add(p)
    members = list(syn2.spine_members)
    for (b, a, c), mv in syn2.Mem.items():
        mk = translate_key(syn2.mkey((b, a, c)), m_new, syn.m)
        v = syn.Mkey.get(mk)
        if v is None:
            raise KeyError(('member', mk))
        if v in ms:
            members.append((a, b, c))
    return syn2, E, members



def spec_instances(spec, syn):
    """return (must_extend list of (sigma, locked), must_not list of (sigma, locked)) for a given instance"""
    m = syn.m
    sig = valid_sigmas(m)
    must, mustnot = [], []
    if spec == 'mexact':
        for s_ in sig:
            (must if s_.count('M') == 1 else mustnot).append((s_, ()))
    elif spec == 'pfg':
        for i in range(m):
            must.append((None, tuple(k for k in range(m) if k != i)))
        mustnot.append((None, tuple(range(m))))
    return must, mustnot


def fix_lists(syn, sigma):
    fa, ft = [], []
    if sigma is None:
        return fa, ft
    for i, st in enumerate(sigma):
        a, b, c = syn.spine_members[i]
        tab, tbc, tac = STATE[st]
        for (p, val) in ((key(a, b), tab), (key(b, c), tbc), (key(a, c), tac)):
            (ft if val else fa).append(p)
    return fa, ft


def run_multi(spec, ms, ncls, cmax, Mmem, seconds, seed, anchors, glue, extend=0):
    shared = Shared()
    syns = [PerSynth(m, ncls, cmax, Mmem, anchors=anchors, glue=glue, shared=shared) for m in ms]
    specs = []
    for syn in syns:
        must, mustnot = spec_instances(spec, syn)
        for (sg, lk) in must:
            syn.add_scenario(sig=sg, locked=lk)
        specs.append((must, mustnot))
    solver = Cadical153(bootstrap_with=shared.cl)
    t0 = time.time()
    it = 0
    while time.time() - t0 < seconds:
        if not solver.solve():
            print('UNSAT after', it, 'cex', flush=True)
            return None
        model = solver.get_model()
        ms_set = set(l for l in model if l > 0)
        it += 1
        bad = None
        for syn, (must, mustnot) in zip(syns, specs):
            E, members = syn.decode(model)
            for (sg, lk) in mustnot:
                fa, ft = fix_lists(syn, sg)
                sing = [syn.anchor_v[k] for k in lk] if lk else []
                lab = ze_solve(syn.n, E, members, singletons=sing, fixed_apart=fa, fixed_together=ft)
                if lab is not None:
                    bad = (syn, lab)
                    break
            if bad:
                break
        if bad is None:
            print('FOUND after', it, flush=True)
            return shared, syns, ms_set
        solver.add_clause(bad[0].cex_clause(bad[1]))
        if it % 100 == 0:
            print('iter', it, round(time.time() - t0, 1), flush=True)
    print('TIMEOUT', it, flush=True)
    return 'timeout'


def evaluate(spec, syn, ms_set, shared):
    class _M:  # emulate model list
        pass
    model = list(ms_set)
    E, members = syn.decode(model)
    ok = check_instance(syn.n, E, members)[0]
    res = dict(m=syn.m, n=syn.n, valid=ok, H=len(members), c=max(sum(1 for mm in members if v in mm) for v in range(syn.n)))
    if spec == 'mexact':
        sols = ze_solve(syn.n, E, members, all_solutions=True, limit=100000)
        res['spine'] = sorted(set(spine_states2(syn, l) for l in sols))
        res['good'] = res['spine'] == sorted(set(''.join(['R'] * a + ['M'] + ['L'] * (syn.m - a - 1)) for a in range(syn.m)))
    elif spec == 'pfg':
        A = syn.anchor_v
        faces = {}
        res['all_locked_feasible'] = ze_solve(syn.n, E, members, singletons=A) is not None
        res['all_but_i'] = [ze_solve(syn.n, E, members, singletons=[x for x in A if x != a]) is not None for a in A]
    return res, E, members


def export_base(syn, ms_set):
    E, members = syn.decode(list(ms_set))
    Enames = set(frozenset((syn.vname[p[0]], syn.vname[p[1]])) for p in E)
    mnames = [(syn.vname[a], syn.vname[b], syn.vname[c]) for (a, b, c) in members]
    return Enames, mnames


def load_base(prefix, m, nbase=5):
    import json
    d = json.load(open('%s_m%d.json' % (prefix, m)))
    vn = [tuple(x) for x in d['vname']]
    Enames = set(frozenset((vn[a], vn[b])) for (a, b) in d['E'])
    mnames = [(vn[a], vn[b], vn[c]) for (a, b, c) in d['members']]
    return dict(E=Enames, members=mnames, nbase=nbase)


def run_addon(prefix, ms, k, cmax, Mmem, seconds, seed, extra_kw=None):
    shared = Shared()
    nb = int(os.environ.get('NBASE', '5'))
    ncls = nb + k
    syns = [PerSynth(m, ncls, cmax, Mmem, anchors=True, glue=True, shared=shared, base=load_base(prefix, m, nb)) for m in ms]
    specs = []
    for syn in syns:
        must, mustnot = spec_instances('pfg', syn)
        for (sg, lk) in must:
            syn.add_scenario(sig=sg, locked=lk)
        specs.append((must, mustnot))
    solver = Cadical153(bootstrap_with=shared.cl)
    t0 = time.time()
    it = 0
    while time.time() - t0 < seconds:
        if not solver.solve():
            print('UNSAT after', it, 'cex', flush=True)
            return None
        model = solver.get_model()
        ms_set = set(l for l in model if l > 0)
        it += 1
        bad = None
        for syn, (must, mustnot) in zip(syns, specs):
            E, members = syn.decode(model)
            for (sg, lk) in mustnot:
                sing = [syn.anchor_v[q] for q in lk]
                lab = ze_solve(syn.n, E, members, singletons=sing)
                if lab is not None:
                    bad = (syn, lab)
                    break
            if bad:
                break
        if bad is None:
            print('FOUND after', it, flush=True)
            return shared, syns, ms_set
        solver.add_clause(bad[0].cex_clause(bad[1]))
        if it % 50 == 0:
            print('iter', it, round(time.time() - t0, 1), flush=True)
    print('TIMEOUT', it, flush=True)
    return 'timeout'


if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'addon':
    prefix = sys.argv[2]
    ms = [int(x) for x in sys.argv[3].split(',')]
    k, cmax, Mmem = map(int, sys.argv[4:7])
    sec = float(sys.argv[7]); seed = int(sys.argv[8]) if len(sys.argv) > 8 else 1
    res = run_addon(prefix, ms, k, cmax, Mmem, sec, seed)
    if res and res != 'timeout':
        shared, syns, ms_set = res
        import json
        for syn in syns:
            r, E, members = evaluate('pfg', syn, ms_set, shared)
            print('CHECK', r, flush=True)
            if os.environ.get('OUT'):
                json.dump(dict(n=syn.n, E=sorted(E), members=members, vname=syn.vname, anchors=syn.anchor_v, m=syn.m, spec='pfg'), open(os.environ['OUT'] + '_m%d.json' % syn.m, 'w'))
        sys.exit(0)


if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'multi':
    spec = sys.argv[2]
    ms = [int(x) for x in sys.argv[3].split(',')]
    ncls, cmax, Mmem = map(int, sys.argv[4:7])
    sec = float(sys.argv[7]); seed = int(sys.argv[8]) if len(sys.argv) > 8 else 1
    anchors = os.environ.get('ANCH', '0') == '1'
    glue = os.environ.get('GLUE', '1') == '1'
    extend = int(os.environ.get('EXTEND', '4'))
    res = run_multi(spec, ms, ncls, cmax, Mmem, sec, seed, anchors, glue)
    if res and res != 'timeout':
        shared, syns, ms_set = res
        import json
        for m2 in sorted(set(ms) | set(range(max(ms) + 1, max(ms) + extend + 1))):
            syn = [x for x in syns if x.m == m2]
            if not syn:
                syn = PerSynth(m2, ncls, cmax, Mmem, anchors=anchors, glue=glue, shared=shared)
            else:
                syn = syn[0]
            r, E, members = evaluate(spec, syn, ms_set, shared)
            print('CHECK', r, flush=True)
            if os.environ.get('OUT'):
                json.dump(dict(n=syn.n, E=sorted(E), members=members, vname=syn.vname, anchors=syn.anchor_v, m=syn.m, spec=spec), open(os.environ['OUT'] + '_m%d.json' % m2, 'w'))
        sys.exit(0)

if __name__ == '__main__':
    spec = sys.argv[1]
    m, ncls, cmax, Mmem = map(int, sys.argv[2:6])
    sec = float(sys.argv[6])
    seed = int(sys.argv[7]) if len(sys.argv) > 7 else 1
    anchors = os.environ.get('ANCH', '0') == '1'
    glue = os.environ.get('GLUE', '1') == '1'
    res = run(spec, m, ncls, cmax, Mmem, sec, seed, anchors, glue)
    if res and res != 'timeout':
        syn, E, members = res
        import json
        if os.environ.get('OUT'):
            json.dump(dict(n=syn.n, E=sorted(E), members=members, vname=syn.vname, anchors=syn.anchor_v, m=syn.m, spec=spec), open(os.environ['OUT'], 'w'))
        print('n', syn.n, 'valid', check_instance(syn.n, E, members))
        print('vertex names (cell,class):', syn.vname)
        print('members', [(syn.vname[a], syn.vname[b], syn.vname[c]) for (a, b, c) in members[syn.m:]])
        if os.environ.get('EXTEND'):
            for m2 in range(syn.m + 1, syn.m + int(os.environ['EXTEND']) + 1):
                try:
                    syn2, E2, mem2 = extend_check(syn, syn.model_set, m2, ncls, cmax, Mmem, glue, anchors)
                except KeyError as e:
                    print('extend m', m2, 'missing key', e); break
                sols2 = ze_solve(syn2.n, E2, mem2, all_solutions=True, limit=100000)
                st = sorted(set(spine_states2(syn2, l) for l in sols2))
                print('EXTENDED m =', m2, 'n =', syn2.n, 'valid', check_instance(syn2.n, E2, mem2)[0], '|H| =', len(mem2), 'c =', max(sum(1 for mm in mem2 if v in mm) for v in range(syn2.n)), '#ZE =', len(sols2), 'spine states', st, flush=True)
        if spec in ('mspine', 'mexact'):
            sols = ze_solve(syn.n, E, members, all_solutions=True, limit=10000)
            print(len(sols), 'ZE solutions; spine state vectors:', sorted(set(spine_states2(syn, l) for l in sols)))



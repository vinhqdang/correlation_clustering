"""Certificate that OH_m has exactly m ZE solutions with pairwise distinct spine words R^a M L^{m-1-a}, for ALL m >= m0.
Uniform encoding of DP states (relative cells) + subset construction + spec automaton.  See w8c_transfer.py for the model."""
import sys, os, itertools, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w8c_lib import *
from w8c_oh import oh_instance

def kind_of(nm): return ('se' if nm[1] % 2 == 0 else 'so') if nm[0] == 's' else 'e'

class Inst:
    def __init__(self, m):
        self.m = m
        n, E, members, names, idx = oh_instance(m)
        self.n, self.E, self.members, self.names, self.idx = n, E, members, names, idx
        self.cell = [(nm[1] // 2 if nm[0] == 's' else nm[1]) for nm in names]
        self.kind = [kind_of(nm) for nm in names]
        self.soft = set()
        for t in members: self.soft |= pairs_of(t)
        self.byc = {c: [v for v in range(n) if self.cell[v] == c] for c in range(m + 1)}
        self.mem_at = {c: [] for c in range(m + 1)}
        for t in members: self.mem_at[max(self.cell[v] for v in t)].append(t)

def rel(v, I, c): return (I.kind[v], I.cell[v] - c)     # relative code of vertex v w.r.t. current cell c

def step_states(I, c, st):
    """st = (tuple of ((kind,relcell), block) over window cells c-2,c-1 in terms of offsets (-2,-1) relative to c), oldblocks) ->
       list of (label, newstate) where newstate is relative to cell c (offsets -1, 0)."""
    win, oldb = st
    # map window entries to actual vertices of cells c-2, c-1
    wd = {}
    for ((kd, off), b) in win:
        cc = c - 1 + off
        cand = [v for v in I.byc.get(cc, []) if I.kind[v] == kd]
        assert len(cand) == 1, (kd, off, cc)
        wd[cand[0]] = b
    blocks = sorted(set(wd.values()))
    nxt = (max(blocks) + 1) if blocks else 0
    new_vs = I.byc[c]
    out = []
    def rec(i, assign, nb):
        if i == len(new_vs):
            yield dict(assign); return
        v = new_vs[i]
        for b in [b for b in blocks if b not in oldb] + [bb for bb in range(nxt, nb)]:
            assign[v] = b
            yield from rec(i + 1, assign, nb)
        assign[v] = nb
        yield from rec(i + 1, assign, nb + 1)
        del assign[v]
    for assign in rec(0, {}, nxt):
        full = dict(wd); full.update(assign)
        ok = True
        for u in new_vs:
            for v in full:
                if v == u or (v in assign and v < u): continue
                p = key(u, v)
                if p in I.soft: continue
                if (full[u] == full[v]) != (p in I.E): ok = False; break
            if not ok: break
        if not ok: continue
        for t in I.mem_at[c]:
            a, b, cc = t
            if (full[a] != full[b]) + (full[b] != full[cc]) + (full[a] == full[cc]) != 1: ok = False; break
        if not ok: continue
        label = ''
        if 1 <= c <= I.m:
            j = c - 1
            sa, sb_, sc = I.idx[('s', 2 * j)], I.idx[('s', 2 * j + 1)], I.idx[('s', 2 * j + 2)]
            tt = (int(full[sa] == full[sb_]), int(full[sb_] == full[sc]), int(full[sa] == full[sc]))
            label = {(1, 0, 0): 'R', (0, 1, 0): 'L', (1, 1, 1): 'M'}.get(tt, '?')
        keep = {v: bb for v, bb in full.items() if I.cell[v] >= c - 1}
        dropped = set(bb for v, bb in full.items() if I.cell[v] < c - 1)
        kb = set(keep.values())
        newold = set(b for b in oldb if b in kb) | (dropped & kb)
        lab = {}; k = 0; entries = []
        for v in sorted(keep, key=lambda x: (I.cell[x], I.kind[x])):
            b = keep[v]
            if b not in lab: lab[b] = k; k += 1
            entries.append(((I.kind[v], I.cell[v] - c), lab[b]))
        out.append((label, (tuple(entries), frozenset(lab[b] for b in newold if b in lab))))
    return out

def spec_next(sp, lab):
    if sp == 'dead' or lab == '': return sp
    if sp == 'pre':
        return {'R': 'pre', 'M': 'post'}.get(lab, 'dead')
    if sp == 'post':
        return 'post' if lab == 'L' else 'dead'

def words_by_spec(m):
    I = Inst(m)
    # subset DP: dict (frozenset(states), spec) -> number of distinct word-prefixes
    cur = {(frozenset([((), frozenset())]), 'pre'): 1}
    dfa_states = set()
    for c in range(0, m + 1):
        nxt = {}
        for (S, sp), mult in cur.items():
            by = {}
            for st in S:
                for (lab, ns) in step_states(I, c, st):
                    by.setdefault(lab, set()).add(ns)
            for lab, nss in by.items():
                key2 = (frozenset(nss), spec_next(sp, lab))
                nxt[key2] = nxt.get(key2, 0) + mult
        cur = nxt
        dfa_states |= set(cur.keys())
    res = {}
    for (S, sp), mult in cur.items():
        res[sp] = res.get(sp, 0) + mult
    return res, len(dfa_states)

if __name__ == '__main__':
    mmax = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    ok = True
    maxd = 0
    for m in range(7, mmax + 1):
        res, d = words_by_spec(m)
        maxd = max(maxd, d)
        good = (res == {'post': m})
        ok &= good
        if m <= 16 or m % 10 == 0: print(m, res, 'DFA states seen', d, 'OK' if good else 'BAD', flush=True)
    print('all OK up to', mmax, ok, 'max #DFA states', maxd)

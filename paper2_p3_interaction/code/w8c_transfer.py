"""Exact transfer-matrix (cell-by-cell DP) for chains of bounded span.  Used to certify, for ALL m, properties of the
one-hot chain OH_m (w8c_oh.py) via the linear-recurrence argument:
   N_m = #ZE partitions of OH_m is the output of a finite automaton with d states, hence N_m satisfies a linear
   recurrence of order <= d; N_m - m therefore satisfies one of order <= d+2 and vanishes for all m >= m0 as soon as it
   vanishes for d+2 consecutive values of m.
Cell of s_k is k//2 (s_{2m} is cell m), cell of e_j is j; every pair at cell distance >= 3 is a hard non-edge (checked)."""
import sys, os, itertools, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w8c_lib import *
from w8c_oh import oh_instance

def build(m):
    n, E, members, names, idx = oh_instance(m)
    cell = [ (nm[1] // 2 if nm[0] == 's' else nm[1]) for nm in names]
    soft = set()
    for t in members: soft |= pairs_of(t)
    for u, v in itertools.combinations(range(n), 2):
        if abs(cell[u] - cell[v]) >= 3:
            assert (u, v) not in soft and (u, v) not in E, ('long pair', names[u], names[v])
    return n, E, members, names, cell, soft

class Layered:
    def __init__(self, m):
        self.m = m
        self.n, self.E, self.members, self.names, self.cell, self.soft = build(m)
        self.byc = {c: [v for v in range(self.n) if self.cell[v] == c] for c in range(m + 1)}
        self.mem_at = {c: [] for c in range(m + 1)}
        for t in self.members:
            self.mem_at[max(self.cell[v] for v in t)].append(t)
    def pair_edge(self, u, v): return key(u, v) in self.E

def step(L, c, states):
    """states: dict state->count where state = (tuple of (vertex, block) for window vertices (cells c-2, c-1), frozenset of old blocks).
    Returns the dict after adding cell c and dropping cell c-2."""
    new_vs = L.byc[c]
    out = {}
    for (win, oldb), cnt in states.items():
        wd = dict(win)           # vertex -> block id
        blocks = sorted(set(wd.values()))
        nxt = (max(blocks) + 1) if blocks else 0
        # enumerate assignments of new vertices
        def rec(i, assign, nb):
            if i == len(new_vs):
                yield dict(assign)
                return
            v = new_vs[i]
            options = [b for b in blocks if b not in oldb] + [bb for bb in range(nxt, nb)]
            # block ids nxt..nb-1 are blocks created by earlier new vertices
            for b in options:
                assign[v] = b
                yield from rec(i + 1, assign, nb)
            assign[v] = nb     # fresh block
            yield from rec(i + 1, assign, nb + 1)
            del assign[v]
        for assign in rec(0, {}, nxt):
            full = dict(wd); full.update(assign)
            ok = True
            # hard pairs with max cell c: togetherness must equal edge-ness unless soft
            for u in new_vs:
                for v, b in full.items():
                    if v == u: continue
                    p = key(u, v)
                    if (v in assign and v < u): continue
                    if p in L.soft: continue
                    tog = full[u] == full[v]
                    if tog != (p in L.E): ok = False; break
                if not ok: break
            if not ok: continue
            for t in L.mem_at[c]:
                a, b, cc = t
                cnt_e = (full[a] != full[b]) + (full[b] != full[cc]) + (full[a] == full[cc])
                if cnt_e != 1: ok = False; break
            if not ok: continue
            # drop cell c-2
            keep = {v: bb for v, bb in full.items() if L.cell[v] >= c - 1}
            dropped_blocks = set(bb for v, bb in full.items() if L.cell[v] < c - 1)
            newold = set(b for b in oldb if b in set(keep.values())) | (dropped_blocks & set(keep.values()))
            # canonical relabel
            lab = {}; k = 0; ccw = []
            for v in sorted(keep):
                b = keep[v]
                if b not in lab: lab[b] = k; k += 1
                ccw.append((v, lab[b]))
            st = (tuple(ccw), frozenset(lab[b] for b in newold if b in lab))
            out[st] = out.get(st, 0) + cnt
    return out

def count(m, upto=None):
    L = Layered(m)
    states = {((), frozenset()): 1}
    sizes = []
    for c in range(0, m + 1):
        states = step(L, c, states)
        sizes.append(len(states))
    return sum(states.values()), sizes

if __name__ == '__main__':
    for m in range(7, 14):
        t = time.time()
        N, sizes = count(m)
        print(m, N, 'max #states', max(sizes), round(time.time() - t, 1))


def step_nfa(L, c, st):
    """one-state version returning list of (label, newstate); label = state letter of P_{c-1} (if c>=1 and c-1<m) else ''."""
    win, oldb = st
    wd = dict(win)
    new_vs = L.byc[c]
    blocks = sorted(set(wd.values()))
    nxt = (max(blocks) + 1) if blocks else 0
    res = []
    def rec(i, assign, nb):
        if i == len(new_vs):
            yield dict(assign); return
        v = new_vs[i]
        options = [b for b in blocks if b not in oldb] + [bb for bb in range(nxt, nb)]
        for b in options:
            assign[v] = b
            yield from rec(i + 1, assign, nb)
        assign[v] = nb
        yield from rec(i + 1, assign, nb + 1)
        del assign[v]
    for assign in rec(0, {}, nxt):
        full = dict(wd); full.update(assign)
        ok = True
        for u in new_vs:
            for v, b in full.items():
                if v == u or (v in assign and v < u): continue
                p = key(u, v)
                if p in L.soft: continue
                if (full[u] == full[v]) != (p in L.E): ok = False; break
            if not ok: break
        if not ok: continue
        for t in L.mem_at[c]:
            a, b, cc = t
            if (full[a] != full[b]) + (full[b] != full[cc]) + (full[a] == full[cc]) != 1: ok = False; break
        if not ok: continue
        label = ''
        if 1 <= c <= L.m:
            j = c - 1
            sa, sb_, sc = (L.names.index(('s', 2 * j)), L.names.index(('s', 2 * j + 1)), L.names.index(('s', 2 * j + 2)))
            tab, tbc, tac = full[sa] == full[sb_], full[sb_] == full[sc], full[sa] == full[sc]
            label = {(1, 0, 0): 'R', (0, 1, 0): 'L', (1, 1, 1): 'M'}.get((int(tab), int(tbc), int(tac)), '?')
        keep = {v: bb for v, bb in full.items() if L.cell[v] >= c - 1}
        dropped_blocks = set(bb for v, bb in full.items() if L.cell[v] < c - 1)
        kb = set(keep.values())
        newold = set(b for b in oldb if b in kb) | (dropped_blocks & kb)
        lab = {}; k = 0; ccw = []
        for v in sorted(keep):
            b = keep[v]
            if b not in lab: lab[b] = k; k += 1
            ccw.append((v, lab[b]))
        res.append((label, (tuple(ccw), frozenset(lab[b] for b in newold if b in lab))))
    return res

def language_count(m):
    """number of distinct accepted spine words and the set of words (subset construction); also #paths."""
    L = Layered(m)
    start = ((), frozenset())
    # NFA layered by cell; build transitions per cell lazily
    layer = {start}
    trans = []   # per cell: dict state -> list of (label,newstate)
    for c in range(0, m + 1):
        tr = {}
        newlayer = set()
        for st in layer:
            lst = step_nfa(L, c, st)
            tr[st] = lst
            for (_, ns) in lst: newlayer.add(ns)
        trans.append(tr)
        layer = newlayer
    # words: DFS over determinized subsets
    words = {}
    def rec(c, subset, word):
        if c == m + 1:
            words[word] = words.get(word, 0) + 1
            return
        by = {}
        for st in subset:
            for (lab, ns) in trans[c][st]:
                by.setdefault(lab, set()).add(ns)
        for lab, nss in by.items():
            rec(c + 1, frozenset(nss), word + lab)
    rec(0, frozenset([start]), '')
    return words

if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'lang':
    for m in range(7, 15):
        W = language_count(m)
        want = set(''.join(['R'] * a + ['M'] + ['L'] * (m - a - 1)) for a in range(m))
        print(m, 'distinct words', len(W), 'all one-hot', set(W) == want)

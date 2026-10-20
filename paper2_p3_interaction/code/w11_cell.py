"""w11_cell.py -- periodic chain search at c<=cmax: cell template (k members, in-cell sharing with multiplicity<=cmax, port vertices between consecutive cells,
cluster labels with offset 0/1, pendants) replicated N times; min certificate (l=0) for several N; report growth."""
import sys, time, random, itertools
from w11_lib import *
from w11_gen import make_inst

def random_cell(rnd, k, cmax, npend_cell, nlab):
    slots = [(i, s) for i in range(k) for s in range(3)]
    rnd.shuffle(slots)
    nports = rnd.choice([1, 1, 2])
    out_slots = []; in_slots = []
    # choose ports: out slots from members, in slots from other members
    cand = slots[:]
    for _ in range(nports):
        o = cand.pop(); i_ = cand.pop()
        out_slots.append(o); in_slots.append(i_)
    groups = []
    shared = set()
    rest = cand
    while rest:
        a = rest.pop()
        if rnd.random() < 0.25: groups.append([a]); continue
        g = [a]
        target = 3 if (cmax >= 3 and rnd.random() < 0.4) else 2
        rnd.shuffle(rest)
        for b in rest[:]:
            if len(g) >= target: break
            if all(b[0] != x[0] and frozenset((b[0], x[0])) not in shared for x in g):
                g.append(b); rest.remove(b)
        for x in g:
            for y in g:
                if x[0] < y[0]: shared.add(frozenset((x[0], y[0])))
        groups.append(g)
    # vertex ids in the cell: groups..., then out ports, in ports (in port of cell j+1 == out port of cell j)
    return dict(k=k, groups=groups, out=out_slots, inn=in_slots, centres=[rnd.randrange(3) for _ in range(k)],
                perm=[rnd.sample(range(3), 3) for _ in range(k)],
                pend=npend_cell, nlab=nlab,
                lab=None)

def build_chain(cell, N, lab_assign, pend_adj, start_adj, end_adj):
    """lab_assign: dict vertex-key -> (label, offset) or None for singleton cluster."""
    vid = {}
    def V(key_):
        if key_ not in vid: vid[key_] = len(vid)
        return vid[key_]
    H = []
    for j in range(N):
        slotv = {}
        for gi, g in enumerate(cell['groups']):
            for x in g: slotv[x] = ('g', j, gi)
        for pi, (o, i_) in enumerate(zip(cell['out'], cell['inn'])):
            slotv[o] = ('p', j + 1, pi)       # out port of cell j = shared vertex p_{j+1,pi}
            slotv[i_] = ('p', j, pi)          # in port of cell j  = p_{j,pi}
        for i in range(cell['k']):
            t = [V(slotv[(i, s)]) for s in range(3)]
            pm = cell['perm'][i]
            t = [t[pm[0]], t[pm[1]], t[pm[2]]]
            H.append(tuple(t))
    nv = len(vid)
    n = nv
    # pendants
    pend_keys = []
    for j in range(N):
        for q in range(cell['pend']):
            pend_keys.append(('w', j, q)); V(('w', j, q))
    for q in range(len(start_adj)): V(('ws', q))
    for q in range(len(end_adj)): V(('we', q))
    n = len(vid)
    # clusters: each vertex key -> cluster id
    cl = [-1 - v for v in range(n)]      # unique negative by default (singleton)
    inv = {v: k_ for k_, v in vid.items()}
    for v in range(n):
        k_ = inv[v]
        la = lab_assign.get(k_[0] + ('' if True else ''), None) if False else None
    # label by (kind, local index): kind 'g' -> ('g',gi); 'p' -> ('p',pi); 'w' -> ('w',q)
    for v in range(n):
        k_ = inv[v]
        if k_[0] == 'g': lk = ('g', k_[2]); jj = k_[1]
        elif k_[0] == 'p': lk = ('p', k_[2]); jj = k_[1]
        elif k_[0] == 'w': lk = ('w', k_[2]); jj = k_[1]
        elif k_[0] == 'ws': lk = ('ws', k_[1]); jj = 0
        else: lk = ('we', k_[1]); jj = N - 1
        la = lab_assign.get(lk)
        if la is None: continue
        lab, off = la
        cl[v] = (lab, jj + off)
    ids = {}
    cl2 = []
    for v in range(n):
        c = cl[v]
        if isinstance(c, int): cl2.append(-1)
        else:
            if c not in ids: ids[c] = len(ids)
            cl2.append(ids[c])
    return nv, H, n - nv, cl2

def random_labels(rnd, cell, nlab, start_n, end_n):
    keys = [('g', gi) for gi in range(len(cell['groups']))] + [('p', pi) for pi in range(len(cell['out']))] + [('w', q) for q in range(cell['pend'])]
    keys += [('ws', q) for q in range(start_n)] + [('we', q) for q in range(end_n)]
    la = {}
    for k_ in keys:
        if rnd.random() < 0.75: la[k_] = (rnd.randrange(nlab), rnd.choice([0, 0, 1]))
    return la

if __name__ == '__main__':
    k = int(sys.argv[1]); cmax = int(sys.argv[2]); samples = int(sys.argv[3]); seed = int(sys.argv[4]); Ns = [int(x) for x in sys.argv[5].split(',')]
    rnd = random.Random(seed)
    t0 = time.time(); found = 0; stat = {}
    for s in range(samples):
        cell = random_cell(rnd, k, cmax, rnd.choice([1, 2]), rnd.choice([2, 3, 4]))
        sn = rnd.choice([1, 2]); en = rnd.choice([1, 2])
        la = random_labels(rnd, cell, cell['nlab'], sn, en)
        sizes = []
        okc = True
        for N in Ns:
            nv, H, npd, cl = build_chain(cell, N, la, None, [0] * sn, [0] * en)
            I = make_inst(nv, H, npd, cl)
            if not (I.valid() and I.cmax() <= cmax and I.pair_maximal()): okc = False; break
            sz, T = min_certificate(I, 0, max_size=30, workers=2)
            sizes.append(sz)
            if sz is None or sz == 0: break
        if not okc: continue
        key_ = tuple(sizes)
        stat[key_] = stat.get(key_, 0) + 1
        if len(sizes) == len(Ns) and None not in sizes and sizes[-1] > sizes[0] and sizes[0] >= 1:
            found += 1
            print('GROWTH', sizes, 'cell', cell, 'labels', la, flush=True)
    print('done', round(time.time() - t0), 's', 'growth cells', found); print(sorted(stat.items(), key=str))

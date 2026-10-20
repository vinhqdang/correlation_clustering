import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w8c_lib import *
from w8c_thick import *

def toy(p=6, with_lock=True):
    B = Builder()
    M1 = [(0, 1), (2, 3), (4, 5)]
    M2 = [(0, 2), (1, 4), (3, 5)]
    assert p == 6
    pads = set()
    blow_up(B, [('a', 'b', 'c')], [], p, ['a', 'b', 'c'])
    for v in 'abc':
        for (i, k) in M1: pads.add(frozenset([(v, i), (v, k)]))
    for v in 'ab':
        for (i, k) in M2: pads.add(frozenset([(v, i), (v, k)]))
    # z: permanent guest forced into A
    for (i, k) in M1:
        B.member(('a', i), 'z', ('a', k))
    for v in 'bc':
        for (i, k) in M1:
            B.member('z', (v, i), (v, k))
    if with_lock:
        for v in 'ab':
            for (i, k) in M2:
                B.member('y', (v, i), (v, k))
        B.member('x', 'y', 'z')
    # hard cliques minus pads
    for v in 'abc':
        for i in range(p):
            for k in range(i + 1, p):
                if frozenset([(v, i), (v, k)]) in pads: continue
                B.hard((v, i), (v, k))
    return B

if __name__ == '__main__':
    B = toy()
    n, E, members = B.finish()
    print('n', n, '|H|', len(members), check_instance(n, E, members))
    nm = B.names; idx = B.idx
    def cl(lab):
        d = {}
        for v, l in enumerate(lab): d.setdefault(l, []).append(nm[v])
        return sorted(d.values(), key=lambda c: (len(c), str(c)))
    def desc(lab):
        # which of a,b,c proto-clusters are together
        def tog(u, w): return lab[idx[(u, 0)]] == lab[idx[(w, 0)]]
        return ('ab' if tog('a','b') else '-') + ('bc' if tog('b','c') else '-') + ('ac' if tog('a','c') else '-'), \
               'z~a' if lab[idx['z']] == lab[idx[('a',0)]] else 'z!~a', \
               'y~z' if lab[idx['y']] == lab[idx['z']] else ('y~x' if lab[idx['y']] == lab[idx['x']] else 'y?'), 'x single' if sum(1 for l in lab if l == lab[idx['x']]) == 1 else 'x not single'
    sols = ze_solve(n, E, members, all_solutions=True, limit=1000)
    print('#ZE', len(sols))
    for lab in sols: print(desc(lab))
    print('locked (x singleton):')
    sols2 = ze_solve(n, E, members, singletons=[idx['x']], all_solutions=True, limit=1000)
    print('#ZE locked', len(sols2))
    for lab in sols2: print(desc(lab))

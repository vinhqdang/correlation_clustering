"""w11_table.py -- verification of the lock table (list induced on member A=(a,b,c) by an outside vertex z with S=N(z)&A, all three z-pairs fixed),
for the 3-state ZE model and the 5-state model; and of the partial guard types (two z-pairs)."""
import itertools
TOG3 = {'L': (0,1,0), 'R': (1,0,0), 'M': (1,1,1)}                      # (ab,bc,ac) togetherness after the edit
TOG5 = dict(TOG3); TOG5['4'] = (0,0,0); TOG5['5'] = (0,0,1)
def table(S, TOG, pairs=((0,1,0),(1,2,1),(0,2,2))):
    """allowed states: no triple {z,x,y} (xy a pair of A) with exactly two edges"""
    out = []
    for s, t in TOG.items():
        ok = True
        for (x, y, k) in pairs:
            e = (1 if x in S else 0) + (1 if y in S else 0) + t[k]
            if e == 2: ok = False
        if ok: out.append(s)
    return ''.join(out)
if __name__ == '__main__':
    names = {frozenset(): 'empty'}
    for S in [(), (0,), (1,), (2,), (0,1), (0,2), (1,2), (0,1,2)]:
        print('S =', ''.join('abc'[i] for i in S) or '-', ' ZE:', table(set(S), TOG3) or 'EMPTY', ' 5-state:', table(set(S), TOG5) or 'EMPTY')
    print('partial guards (two z-pairs xy):')
    for pair, nm in (((0,1,0),'ab'),((1,2,1),'bc'),((0,2,2),'ac')):
        for sx, sy in ((1,0),(0,1),(1,1),(0,0)):
            S = set()
            if sx: S.add(pair[0])
            if sy: S.add(pair[1])
            print(' pair', nm, 'z-adjacency', (sx, sy), 'list', table(S, TOG3, pairs=(pair,)))

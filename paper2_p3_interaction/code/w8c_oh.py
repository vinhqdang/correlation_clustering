"""The one-hot selection chain OH_m (found by SAT/CEGIS, see notes/wave8_hardness.md Sec. 3).
Vertices s_0..s_{2m} (spine path), e_0..e_m.  Names: ('s',k), ('e',j).  Members listed as (a,b,c) with centre b.
Interior (2<=j<=m-4): G_j=(e_j,s_{2j+1},e_{j+2});  (2<=j<=m-3): F_j=(e_j,s_{2j+2},e_{j+1});  spine P_j=(s_{2j},s_{2j+1},s_{2j+2}).
Caps are explicit below.  Graph edges = member edges (ab,bc) plus the 'hard' list; every other pair is a hard non-edge."""
import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w8c_lib import *

def oh_instance(m):
    assert m >= 7
    S = lambda k: ('s', k)
    Ev = lambda j: ('e', j)
    mem = []
    for j in range(m):
        mem.append((S(2 * j), S(2 * j + 1), S(2 * j + 2)))                 # P_j
    # left cap
    mem += [(Ev(0), S(2), Ev(1)), (S(1), Ev(1), S(3)), (Ev(0), S(4), Ev(2))]
    for j in range(2, m - 3):
        mem.append((Ev(j), S(2 * j + 1), Ev(j + 2)))                        # G_j  (j<=m-4)
    for j in range(2, m - 2):
        mem.append((Ev(j), S(2 * j + 2), Ev(j + 1)))                        # F_j  (j<=m-3)
    # right cap
    mem += [(S(2 * m - 5), Ev(m - 3), S(2 * m - 1)), (Ev(m - 2), S(2 * m - 2), Ev(m)), (S(2 * m - 4), Ev(m - 1), Ev(m)),
            (S(2 * m - 3), Ev(m - 1), S(2 * m)), (Ev(m - 2), Ev(m - 1), S(2 * m - 1)), (S(2 * m - 3), Ev(m), S(2 * m - 1))]
    hard = [(S(0), Ev(1)), (Ev(0), S(3)), (S(2 * m - 3), Ev(m - 2)), (S(2 * m), Ev(m))]
    names = [('s', k) for k in range(2 * m + 1)] + [('e', j) for j in range(m + 1)]
    idx = {v: i for i, v in enumerate(names)}
    members = [tuple(idx[v] for v in t) for t in mem]
    E = set()
    for (a, b, c) in members:
        E.add(key(a, b)); E.add(key(b, c))
    for (u, v) in hard:
        E.add(key(idx[u], idx[v]))
    return len(names), E, members, names, idx

def spine_state_vector(m, idx, lab):
    out = []
    for j in range(m):
        a, b, c = idx[('s', 2 * j)], idx[('s', 2 * j + 1)], idx[('s', 2 * j + 2)]
        out.append({(1, 0, 0): 'R', (0, 1, 0): 'L', (1, 1, 1): 'M'}.get((int(lab[a] == lab[b]), int(lab[b] == lab[c]), int(lab[a] == lab[c])), '?'))
    return ''.join(out)

if __name__ == '__main__':
    for m in range(7, 17):
        n, E, members, names, idx = oh_instance(m)
        ok = check_instance(n, E, members)
        sols = ze_solve(n, E, members, all_solutions=True, limit=100000)
        sv = sorted(spine_state_vector(m, idx, l) for l in sols)
        good = sv == sorted(''.join(['R'] * a + ['M'] + ['L'] * (m - a - 1)) for a in range(m))
        c = max(sum(1 for t in members if v in t) for v in range(n))
        print('m', m, 'n', n, '|H|', len(members), 'valid', ok[0], 'c', c, '#ZE', len(sols), 'one-hot', good)

"""Chain of 'filler' gadgets (found by w8c_det.py, MODE=forbid FMODE=3 q=3 j=1) on a plain P3 spine.
Each interior even vertex s_{2j} gets twins e8 (of s_2j), e9 (of s_{2j+1}), e11 (of s_{2j-1}) and members
M1=(e9,s2j,e11) M2=(s_{2j-1},e8,s_{2j+1}) M3=(s_{2j+2},e9,e8) M4=(s_{2j-2},e11,e8)."""
import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w8c_lib import *

def build_chain(m, positions):
    """spine with m members (vertices 0..2m); filler gadgets at the even vertices s_{2j}, j in positions (1<=j<=m-1)."""
    n = 2 * m + 1
    E = set(key(i, i + 1) for i in range(2 * m))
    members = [(2 * i, 2 * i + 1, 2 * i + 2) for i in range(m)]
    names = {}
    for j in positions:
        e8, e9, e11 = n, n + 1, n + 2
        n += 3
        names[j] = (e8, e9, e11)
        s = lambda k: 2 * j - 2 + k      # window index -> global spine index
        E |= {key(e11, s(0)), key(e11, s(1)), key(e11, s(2)), key(e8, s(1)), key(e8, s(2)), key(e8, s(3)),
              key(e9, s(2)), key(e9, s(3)), key(e9, s(4)), key(e8, e9), key(e8, e11)}
        members += [(e9, s(2), e11), (s(1), e8, s(3)), (s(4), e9, e8), (s(0), e11, e8)]
    return n, E, members, names

def spine_states(m, lab):
    out = []
    for i in range(m):
        a, b, c = 2 * i, 2 * i + 1, 2 * i + 2
        tab, tbc, tac = lab[a] == lab[b], lab[b] == lab[c], lab[a] == lab[c]
        out.append({(1, 0, 0): 'R', (0, 1, 0): 'L', (1, 1, 1): 'M'}.get((int(tab), int(tbc), int(tac)), '?'))
    return ''.join(out)

if __name__ == '__main__':
    m = int(sys.argv[1]); pos = [int(x) for x in sys.argv[2].split(',')]
    n, E, members, names = build_chain(m, pos)
    print('valid', check_instance(n, E, members), 'n', n, '|H|', len(members))
    sols = ze_solve(n, E, members, all_solutions=True, limit=100000)
    sig = sorted(set(spine_states(m, l) for l in sols))
    print(len(sols), 'ZE solutions; spine state vectors:', sig)

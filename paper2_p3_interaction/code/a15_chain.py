import sys
from w8c_oh import *
for m in [int(x) for x in sys.argv[1:]]:
    n, E, members, names, idx = oh_instance(m)
    sols = ze_solve(n, E, members, all_solutions=True, limit=1000)
    sv = sorted(spine_state_vector(m, idx, l) for l in sols)
    print('m', m, '#ZE', len(sols), 'one-hot', sv == sorted(''.join(['R'] * a + ['M'] + ['L'] * (m - a - 1)) for a in range(m)), flush=True)

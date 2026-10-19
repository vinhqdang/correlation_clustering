"""control for Theorem P: use threshold l-1 in the pruning (should FAIL somewhere if the test is sensitive)."""
import sys, random
import a7_lib as L, a7_thmP as T
import networkx as nx
rng = random.Random(int(sys.argv[1])); fails = 0; tot = 0
for it in range(int(sys.argv[2])):
    n = rng.choice([6,7,8]); A,H = L.gen_instance(n, rng)
    if not H or not L.is_pair_maximal(A,H): continue
    opt = L.opt_subset(A)
    for l in range(1,4):
        st = T.setup(A,H,l-1)   # prune with threshold l-1 but budget l
        st['Z'] = 3*st['D']+3*l+1
        m_all, m_conn = T.brute_m(A, st)
        sumK = sum(T.f_true(A,K) for K in st['comps'] if any(v in st['BIG'] for v in K))
        base = int(A.sum()//2)+sumK
        tot += 1
        if (base+m_conn <= len(H)+l) != (opt <= len(H)+l): fails += 1
print("control: instances-l pairs", tot, "criterion failures with threshold l-1:", fails)

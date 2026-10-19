"""l = 0 dynamic program on the real PACE instances, using the local-search packing of
exp_packings.py (same seeds).  Records the elimination width of the CONSTRAINT graph
(a subgraph of I(H): an edge only where some non-trivial triple constraint links two
P3s), runs the DP when that width is <= WMAX, and checks DP answer == (OPT == |H|).
    python3 exp_dp.py   -> results_dp.csv
"""
import csv, os, sys, time
import numpy as np
import pandas as pd
from multiprocessing import Pool
from exp_common import *
from exp_solver import build_csp, constraint_graph, dp_l0, unpacked_p3

HERE = os.path.dirname(os.path.abspath(__file__))
WMAX = 16
HMAX = 3000

def job(args):
    fn, opt = args
    n, edges = read_pace(os.path.join(ROOT, 'data/raw/_pace21/exact', fn))
    A = adjacency(n, edges); indptr, indices = csr(A)
    Hg = greedy_pack(A, indptr, indices, 0, n + 5)
    iters = min(300000, 100 * len(Hg) + 20000)
    H, _ = local_search(A, Hg, iters, 1, 4)
    row = dict(instance=fn, n=n, p=len(H), opt=opt, opt_eq_H=int(opt == len(H)))
    if len(H) > HMAX:
        return row
    t0 = time.time()
    cons = build_csp(A, H)
    adj = constraint_graph(len(H), cons)
    row['ncons'] = len(cons); row['csp_sec'] = round(time.time() - t0, 2)
    row['cg_maxdeg'] = max(len(a) for a in adj)
    row['cg_edges'] = sum(len(a) for a in adj) // 2
    w1, _ = elimination(adj, 150, 'mindeg')
    w2, _ = elimination(adj, 150, 'minfill') if len(H) <= 1500 else (None, None)
    ws = [w for w in (w1, w2) if w is not None]
    row['cg_tw_ub'] = min(ws) if ws else None
    if ws and min(ws) <= WMAX:
        t1 = time.time()
        ans, info = dp_l0(A, H, max_width=WMAX, check_unpacked=False)
        row['dp_answer'] = ans; row['dp_sec'] = round(time.time() - t1, 2)
        row['dp_correct'] = int(ans == bool(row['opt_eq_H']))
    return row

if __name__ == '__main__':
    d = pd.read_csv(os.path.join(ROOT, 'data/pace2021_exact_kapoce_bounds.csv'))
    d = d[d.solved_1h == 1]
    jobs = [(r.instance, int(r.opt_or_empty)) for r in d.itertuples()]
    with Pool(4) as p:
        rows = p.map(job, jobs, chunksize=1)
    keys = []
    for r in rows:
        for k in r:
            if k not in keys: keys.append(k)
    with open(os.path.join(HERE, 'results_dp.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)
    print('done')

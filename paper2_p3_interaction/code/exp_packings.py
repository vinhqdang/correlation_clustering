"""P3 packings (greedy + local search) and interaction statistics on real graphs.

  python3 exp_packings.py pace    -> results_pace.csv   (173 solved PACE 2021 exact instances)
  python3 exp_packings.py snap    -> results_snap.csv   (small SNAP graphs)
Seeds: greedy seed 0, local search seed 1.  Deterministic given numba's RNG.
"""
import csv, sys, os, time
import numpy as np
from multiprocessing import Pool
from exp_common import *
from exp_solver import dp_l0

HERE = os.path.dirname(os.path.abspath(__file__))
SNAP = ['email-Eu-core', 'soc-sign-bitcoinalpha', 'facebook_combined', 'soc-sign-bitcoinotc', 'ca-GrQc', 'wiki-Vote', 'ca-HepTh']
SNAPFILE = {'email-Eu-core': ('email-Eu-core.txt.gz', -1), 'soc-sign-bitcoinalpha': ('soc-sign-bitcoinalpha.csv.gz', 2),
            'soc-sign-bitcoinotc': ('soc-sign-bitcoinotc.csv.gz', 2), 'facebook_combined': ('facebook_combined.txt.gz', -1),
            'ca-GrQc': ('ca-GrQc.txt.gz', -1), 'wiki-Vote': ('wiki-Vote.txt.gz', -1), 'ca-HepTh': ('ca-HepTh.txt.gz', -1)}

def analyse(name, n, edges, opt, do_dp, ls_iters_cap=300000):
    t0 = time.time()
    A = adjacency(n, edges)
    m = int(A.sum() // 2)
    indptr, indices = csr(A)
    big = n + 5
    row = dict(instance=name, n=n, m=m, density=round(m / (n * (n - 1) / 2), 5),
               maxdeg_G=int(A.sum(1).max()), opt=opt)
    Hg = greedy_pack(A, indptr, indices, 0, big)
    row['H_greedy'] = len(Hg)
    for cap in (1, 2, 3):
        row[f'H_cap{cap}'] = len(greedy_pack(A, indptr, indices, 0, cap))
    iters = min(ls_iters_cap, 100 * len(Hg) + 20000)
    t1 = time.time()
    H, gains = local_search(A, Hg, iters, 1, 4)
    row['ls_iters'] = iters; row['ls_sec'] = round(time.time() - t1, 1)
    row['H_ls'] = len(H)
    if n <= 700:
        verify_packing(A, H); row['maximal'] = is_maximal(A, H)
    for tag, HH in (('g', Hg), ('ls', H)):
        s = interaction_stats(HH, n)
        for k in ('c', 'comps', 'maxdeg', 'avgdeg', 'tw_lb', 'tw_ub', 'tw_ub_mindeg', 'tw_ub_minfill'):
            row[f'{k}_{tag}'] = s[k]
    if opt is not None:
        row['ratio_g'] = round(len(Hg) / opt, 4); row['ratio_ls'] = round(len(H) / opt, 4)
        row['opt_eq_H'] = int(opt == len(H))
    if do_dp:
        t2 = time.time()
        if row['H_ls'] <= 4000:
            ans, info = dp_l0(A, H, max_width=14, check_unpacked=False)
            row['dp_answer'] = ans; row['dp_width'] = info['width']; row['dp_ncons'] = info['ncons']
            row['dp_sec'] = round(time.time() - t2, 2)
    row['total_sec'] = round(time.time() - t0, 1)
    return row

def pace_job(args):
    fn, opt = args
    n, edges = read_pace(os.path.join(ROOT, 'data/raw/_pace21/exact', fn))
    return analyse(fn, n, edges, opt, True)

def snap_job(name):
    fn, sc = SNAPFILE[name]
    n, edges = read_snap(os.path.join(ROOT, 'data/raw', fn), sc)
    return analyse(name, n, edges, None, False, 200000)

def write(rows, out):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)

if __name__ == '__main__':
    which = sys.argv[1]
    if which == 'pace':
        import pandas as pd
        d = pd.read_csv(os.path.join(ROOT, 'data/pace2021_exact_kapoce_bounds.csv'))
        d = d[d.solved_1h == 1]
        jobs = [(r.instance, int(r.opt_or_empty)) for r in d.itertuples()]
        if len(sys.argv) > 2:
            jobs = jobs[:int(sys.argv[2])]
        with Pool(4) as p:
            rows = p.map(pace_job, jobs, chunksize=1)
        write(rows, os.path.join(HERE, 'results_pace.csv'))
    else:
        with Pool(3) as p:
            rows = p.map(snap_job, SNAP, chunksize=1)
        write(rows, os.path.join(HERE, 'results_snap.csv'))
    print('done', len(rows))

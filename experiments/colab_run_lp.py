"""Relaxation strength and baselines for the lower bounds (Section 5).

Modes (the tag's second character onward, e.g. tag "ltri"):
  pack  sparse star packing improved by ruin-and-recreate local search for T
        seconds: the scalable form of KaPoCE's root star-packing bound; written
        as a certificate and checked like every other bound;
  tri   metric LP restricted to P (triangle rows only), cutting planes until
        no violated row or T seconds; LP value (floating point, not certified);
  star  the same with star rows;
  full  one block covering the whole graph, triangle + star + local subgraph
        rows, cold start, until convergence or T; certificate checked;
  warm  as full, after the star packing warm start (the CertiFlip bound with no
        time limit); certificate checked.

usage: colab_run_lp.py T SEED TAG GRAPH..."""
import json, os, platform, shutil, subprocess, sys, time
CC = os.environ.get('CC_ROOT', '/content/cc')
sys.path.insert(0, CC); sys.path.insert(0, os.path.join(CC, 'experiments'))
OUT = os.environ.get('RESULTS_DIR', '/content/results')
import numpy as np


def cpu_model():
    try:
        for line in subprocess.run(['lscpu'], capture_output=True, text=True).stdout.splitlines():
            if line.startswith('Model name'):
                return line.split(':', 1)[1].strip()
    except OSError:
        pass
    return 'unknown'


def run(name, mode, T, seed):
    import datasets as D
    from ccbench.support import build_support
    from ccbench.blockdual import BlockDualBound, add_star_packing
    from ccbench.certiflip import block_bound, _packing_rate
    from ccbench.dual import star_packing_ls
    from ccbench.lp import SparseLP
    g = D.load(name)
    t0 = time.time()
    sup = build_support(g)
    out = {'graph': name, 'n': g.n, 'm': g.m, 'pairs': int(sup.npairs), 'mode': mode, 'T': T,
           'seed': seed, 'support_time': time.time() - t0}
    bd = None
    t1 = time.time()
    if mode == 'pack':
        bd = BlockDualBound(g, sup)
        _, rate = _packing_rate(g, bd)
        v, rp, rpairs, rk = star_packing_ls(g, sup, pgraph=(bd.ptr, bd.idx, bd.pid),
                                            iters=int(max(1000, rate * T)), seed=seed)
        add_star_packing(bd, rp, rpairs, rk)
        out.update({'value': float(bd.bound()), 'stars': int(len(rk))})
    elif mode in ('tri', 'star'):
        lp = SparseLP(g, sup)
        if mode == 'tri':
            lp.solve(max_rounds=100000, time_limit=T)
        else:
            lp.solve_with_stars(max_rounds=100000, time_limit=T)
        out.update({'value': float(lp.value), 'converged': bool(lp.converged),
                    'rounds': int(lp.rounds), 'rows': int(lp.nrows)})
    elif mode == 'full':
        bd = BlockDualBound(g, sup)
        bd.solve_block(np.arange(g.n), time_limit=T, max_rounds=100000, subgraphs=True)
        out.update({'value': float(bd.bound()),
                    'converged': bool(getattr(bd, 'last_block', {}).get('converged', False))})
    elif mode == 'warm':
        bd = block_bound(g, sup, time_limit=T, block_size=max(2500, g.n), seed=seed)
        out.update({'value': float(bd.bound()),
                    'converged': bool(getattr(bd, 'last_block', {}).get('converged', False))})
    else:
        raise ValueError(mode)
    out['time'] = time.time() - t1
    if bd is not None:
        import check_certificate as C
        import instance_meta as M
        safe = name.replace('/', '_')
        cert = os.path.join(OUT, 'certs', f'l{mode}_{safe}_{seed}.npz')
        os.makedirs(os.path.dirname(cert), exist_ok=True)
        np.savez_compressed(cert, **bd.certificate())
        M.stamp(cert, name, g)
        try:
            rep = C.check_file(M.raw_file(name)[1], cert)
            out.update({'check': 'ok', 'certified': int(rep['certified']),
                        'lb_exact': rep['lb_exact'], 'identity': rep['identity']})
        except ValueError as exc:
            out['check'] = f'rejected: {exc}'
    return out


def main(graphs, T, seed, tag):
    os.makedirs(OUT, exist_ok=True)
    mode = tag.split('+')[0][1:]
    try:
        commit = open(os.path.join(CC, 'COMMIT')).read().strip()
    except OSError:
        commit = 'unknown'
    for name in graphs:
        out = run(name, mode, T, seed)
        out.update({'tag': tag, 'commit': commit, 'cpu': cpu_model(),
                    'python': platform.python_version()})
        safe = name.replace('/', '_')
        with open(os.path.join(OUT, f'{tag}_{safe}_{seed}.json'), 'w') as fh:
            json.dump(out, fh)
        print(name, out, flush=True)


if __name__ == '__main__':
    main(sys.argv[4:], float(sys.argv[1]), int(sys.argv[2]), sys.argv[3])

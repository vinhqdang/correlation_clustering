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
        time limit); certificate checked;
  plp   a long packing followed by block LPs, as one procedure: the packing of
        mode pack with T/2 seconds of iterations, then METIS block sweeps from
        it for T/2 seconds (Theorem 5(b): never below the packing); the value
        after the packing is recorded as well; certificate checked;
  kroot the root star-packing and P3-packing bounds of the KaPoCE
        branch-and-bound (its `lbounds` program, branch `experiments`, built
        by experiments/kapoce/build_lbounds.sh), limit T seconds; integers
        computed by KaPoCE, not certified here.

usage: colab_run_lp.py T SEED TAG GRAPH..."""
import json, os, platform, shutil, subprocess, sys, time
CC = os.environ.get('CC_ROOT', '/content/cc')
sys.path.insert(0, CC); sys.path.insert(0, os.path.join(CC, 'experiments'))
OUT = os.environ.get('RESULTS_DIR', '/content/results')
import numpy as np


def ensure_pace():
    """PACE 2021 instances on this machine (the certificate jobs may run on a
    machine that never ran a head-to-head job)."""
    d = os.path.join(CC, 'data', 'raw', 'pace')
    if all(os.path.isdir(os.path.join(d, t)) and len(os.listdir(os.path.join(d, t))) >= 200
           for t in ('exact', 'heur')):
        return
    sh = lambda c: subprocess.run(c, shell=True, capture_output=True, text=True)
    sh('cd /tmp && rm -rf pace21 && git clone -q --depth 1 '
       'https://github.com/PACE-challenge/Cluster-Editing-PACE-2021-instances pace21')
    for track in ('exact', 'heur'):
        t = os.path.join(d, track)
        os.makedirs(t, exist_ok=True)
        sh(f"cd /tmp/pace21 && for f in $(find . -name '{track}*.gr*'); do cp $f {t}/; done; "
           f"cd {t} && for f in *.gz; do [ -e \"$f\" ] && gunzip -f \"$f\"; done; "
           f"for f in *.xz; do [ -e \"$f\" ] && unxz -f \"$f\"; done")


def cpu_model():
    try:
        for line in subprocess.run(['lscpu'], capture_output=True, text=True).stdout.splitlines():
            if line.startswith('Model name'):
                return line.split(':', 1)[1].strip()
    except OSError:
        pass
    return 'unknown'


def kapoce_lb():
    """The lbounds program, built once per machine (the build is incremental)."""
    import fcntl
    d = os.environ.get('KAPOCE_LB_DIR', '/content/kapoce_lb')
    with open(d.rstrip('/') + '.lock', 'w') as lock:  # one build at a time
        fcntl.flock(lock, fcntl.LOCK_EX)
        r = subprocess.run(['bash', os.path.join(CC, 'experiments', 'kapoce', 'build_lbounds.sh'),
                            d], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError('lbounds build failed: ' + r.stderr[-500:])
    return r.stdout.strip().splitlines()[-1]


def kroot(g, name, T):
    import tempfile
    from baselines import write_pace
    out = {'graph': name, 'n': g.n, 'm': g.m, 'mode': 'kroot', 'T': T, 'kapoce': '63079a9'}
    binary = kapoce_lb()
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, 'g.gr')
        write_pace(g, f)
        t0 = time.time()
        with open(f) as fin:
            try:
                p = subprocess.run([binary], stdin=fin, capture_output=True, text=True,
                                   timeout=T)
                text = p.stdout
                out['status'] = 'ok' if p.returncode == 0 else f'exit {p.returncode}'
            except subprocess.TimeoutExpired as exc:
                text = exc.stdout or ''
                text = text.decode() if isinstance(text, bytes) else text
                out['status'] = 'timeout'
    out['time'] = time.time() - t0
    # "p3 LB SECONDS" then "star LB SECONDS"; the P3 line survives a timeout
    for line in text.splitlines():
        p = line.split()
        if len(p) == 3 and p[0] in ('p3', 'star'):
            out[p[0]], out[p[0] + '_time'] = int(p[1]), float(p[2])
    out['value'] = out.get('star')
    return out


def run(name, mode, T, seed):
    import datasets as D
    from ccbench.support import build_support
    from ccbench.blockdual import BlockDualBound, add_star_packing
    from ccbench.certiflip import block_bound, _packing_rate
    from ccbench.dual import star_packing_ls
    from ccbench.lp import SparseLP
    g = D.load(name)
    if mode == 'kroot':
        return kroot(g, name, T)
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
    elif mode == 'plp':
        bd = BlockDualBound(g, sup)
        _, rate = _packing_rate(g, bd)
        v, rp, rpairs, rk = star_packing_ls(g, sup, pgraph=(bd.ptr, bd.idx, bd.pid),
                                            iters=int(max(1000, rate * T / 2)), seed=seed)
        add_star_packing(bd, rp, rpairs, rk)
        out.update({'pack_value': float(bd.bound()), 'pack_time': time.time() - t1,
                    'stars': int(len(rk))})
        block_bound(g, sup, time_limit=T / 2, seed=seed, bd=bd, packing_fraction=0.0)
        out.update({'value': float(bd.bound())})
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
        r = subprocess.run(['git', '-C', CC, 'rev-parse', '--short', 'HEAD'],
                           capture_output=True, text=True)
        commit = r.stdout.strip() or 'unknown'
    if any(n.startswith('pace-') for n in graphs):
        ensure_pace()
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

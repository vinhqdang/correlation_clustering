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
        computed by KaPoCE, not certified here;
  kstar the same star packing of KaPoCE (limit T seconds), with its stars
        written out and installed as a dual solution on P like our own
        packing; certificate checked;
  sstar our sparse implementation of that star packing heuristic
        (experiments/sstar/sstar.cc) for T seconds; the stars are written
        directly as certificate rows, without building the support, so it also
        runs where the support does not fit into memory; certificate checked;
  sslp  sstar for T/2 seconds, then METIS block sweeps from its packing for
        T/2 seconds (as plp); certificate checked.

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


def kapoce_stars(g, sup, T):
    """The root star packing of KaPoCE as (rp, rpairs, rk) on the pair ids of
    the support, as returned by our own packing, and KaPoCE's bound and time."""
    import tempfile
    from baselines import write_pace
    binary = kapoce_lb()
    with tempfile.TemporaryDirectory() as d:
        f, sf = os.path.join(d, 'g.gr'), os.path.join(d, 'stars.txt')
        write_pace(g, f)
        with open(f) as fin:
            p = subprocess.run([binary, sf], stdin=fin, capture_output=True, text=True,
                               timeout=T)
        line = [l.split() for l in p.stdout.splitlines() if l.startswith('star ')]
        if p.returncode != 0 or not line:
            raise RuntimeError(f'lbounds failed (exit {p.returncode}): {p.stderr[-300:]}')
        value, secs = int(line[-1][1]), float(line[-1][2])
        stars = [list(map(int, l.split())) for l in open(sf) if l.strip()]
    n = g.n
    key = sup.pu.astype(np.int64) * n + sup.pv.astype(np.int64)
    order = np.argsort(key)
    skey = key[order]

    def pid(a, b):
        k = min(a, b) * n + max(a, b)
        i = np.searchsorted(skey, k)
        if i == len(skey) or skey[i] != k:
            raise ValueError(f'pair ({a}, {b}) of a KaPoCE star is not in P')
        return int(order[i])
    rp, rpairs, rk = [0], [], []
    for st in stars:
        c, leaves = st[0], st[1:]
        # centre pairs first, then the pairs between leaves (as in our packing)
        rpairs += [pid(c, l) for l in leaves]
        rpairs += [pid(leaves[i], leaves[j]) for i in range(len(leaves))
                   for j in range(i + 1, len(leaves))]
        rp.append(len(rpairs))
        rk.append(len(leaves))
    return (np.array(rp, dtype=np.int64), np.array(rpairs, dtype=np.int64),
            np.array(rk, dtype=np.int64), value, secs, len(stars))


def sstar_bin():
    """experiments/sstar/sstar, compiled once per machine."""
    import fcntl
    src = os.path.join(CC, 'experiments', 'sstar', 'sstar.cc')
    b = os.environ.get('SSTAR_BIN', os.path.join(os.path.dirname(OUT.rstrip('/')) or '.', 'sstar_bin'))
    with open(b + '.lock', 'w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if not os.path.exists(b) or os.path.getmtime(b) < os.path.getmtime(src):
            r = subprocess.run(['g++', '-O2', '-std=c++17', '-o', b + '.tmp', src],
                               capture_output=True, text=True)
            if r.returncode != 0:
                raise RuntimeError('sstar build failed: ' + r.stderr[-500:])
            os.replace(b + '.tmp', b)
    return b


def sstar_stars(g, T, seed, min_time):
    """Run sstar on g; returns (stars as lists centre-first, value, seconds, init)."""
    import tempfile
    from baselines import write_pace
    binary = sstar_bin()
    with tempfile.TemporaryDirectory() as d:
        f, sf = os.path.join(d, 'g.gr'), os.path.join(d, 'stars.txt')
        write_pace(g, f)
        with open(f) as fin:
            p = subprocess.run([binary, sf, str(T), str(seed), str(min_time)], stdin=fin,
                               capture_output=True, text=True, timeout=T + 3600)
        res = {l.split()[0]: l.split()[1:] for l in p.stdout.splitlines() if l.strip()}
        if p.returncode != 0 or 'star' not in res:
            raise RuntimeError(f'sstar failed (exit {p.returncode}): {p.stderr[-300:]}')
        stars = [list(map(int, l.split())) for l in open(sf) if l.strip()]
    return (stars, int(res['star'][0]), float(res['star'][1]), int(res['init'][0]),
            [int(x) for x in res.get('rounds', [0, 0])])


def star_rows(stars):
    """Rows of a star packing as vertex pairs: centre pairs (+1), then leaf
    pairs (-1); b = (k - 1) - k(k - 1)/2 and y = 1 (as add_star_packing)."""
    us, vs, vals, ptr, bs = [], [], [], [0], []
    for st in stars:
        c, L = st[0], np.asarray(st[1:], dtype=np.int64)
        k = len(L)
        iu, ju = np.triu_indices(k, 1)
        us.append(np.concatenate([np.full(k, c, dtype=np.int64), L[iu]]))
        vs.append(np.concatenate([L, L[ju]]))
        vals.append(np.concatenate([np.ones(k), -np.ones(len(iu))]))
        ptr.append(ptr[-1] + k + len(iu))
        bs.append((k - 1) - k * (k - 1) / 2)
    u, v = np.concatenate(us), np.concatenate(vs)
    return np.minimum(u, v), np.maximum(u, v), np.concatenate(vals), np.array(ptr, dtype=np.int64), \
        np.array(bs, dtype=np.float64)


def run(name, mode, T, seed, tag=None):
    import datasets as D
    from ccbench.support import build_support
    from ccbench.blockdual import BlockDualBound, add_star_packing
    from ccbench.certiflip import block_bound, _packing_rate
    from ccbench.dual import star_packing_ls
    from ccbench.lp import SparseLP
    g = D.load(name)
    if mode == 'kroot':
        return kroot(g, name, T)
    if mode == 'sstar':
        return run_sstar(g, name, T, seed, tag)
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
    elif mode == 'kstar':
        bd = BlockDualBound(g, sup)
        rp, rpairs, rk, value, secs, ns = kapoce_stars(g, sup, T)
        add_star_packing(bd, rp, rpairs, rk)
        out.update({'value': float(bd.bound()), 'kapoce_star': value, 'kapoce_time': secs,
                    'stars': ns, 'kapoce': '63079a9'})
    elif mode == 'sslp':
        bd = BlockDualBound(g, sup)
        stars, value, secs, init, rounds = sstar_stars(g, T / 2, seed, T / 2)
        u, v, vals, ptr, b = star_rows(stars)
        key = sup.pu.astype(np.int64) * g.n + sup.pv.astype(np.int64)
        order = np.argsort(key)
        pos = np.searchsorted(key[order], u * g.n + v)
        if (pos >= len(key)).any() or (key[order][np.minimum(pos, len(key) - 1)] != u * g.n + v).any():
            raise ValueError('a pair of a star is not in P')
        rpairs = order[pos].astype(np.int64)
        rk = np.array([len(st) - 1 for st in stars], dtype=np.int64)
        add_star_packing(bd, ptr, rpairs, rk)
        out.update({'pack_value': float(bd.bound()), 'pack_time': time.time() - t1,
                    'sstar_value': value, 'sstar_time': secs, 'sstar_init': init,
                    'rounds': rounds, 'stars': len(stars)})
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
        # named after the tag (e.g. lpack+eq), so that variants of a mode do not
        # overwrite each other's certificates
        cert = os.path.join(OUT, 'certs', f'{tag or "l" + mode}_{safe}_{seed}.npz')
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


def write_check(cert_d, name, g, tag, seed, out):
    import check_certificate as C
    import instance_meta as M
    safe = name.replace('/', '_')
    cert = os.path.join(OUT, 'certs', f'{tag}_{safe}_{seed}.npz')
    os.makedirs(os.path.dirname(cert), exist_ok=True)
    np.savez_compressed(cert, **cert_d)
    M.stamp(cert, name, g)
    try:
        rep = C.check_file(M.raw_file(name)[1], cert)
        out.update({'check': 'ok', 'certified': int(rep['certified']),
                    'lb_exact': rep['lb_exact'], 'identity': rep['identity']})
    except ValueError as exc:
        out['check'] = f'rejected: {exc}'
    return out


def run_sstar(g, name, T, seed, tag):
    """The sparse star packing alone, its stars written directly as rows: the
    support is never built."""
    out = {'graph': name, 'n': g.n, 'm': g.m, 'mode': 'sstar', 'T': T, 'seed': seed}
    t1 = time.time()
    stars, value, secs, init, rounds = sstar_stars(g, T, seed, T)
    u, v, vals, ptr, b = star_rows(stars)
    out.update({'value': float(value), 'sstar_time': secs, 'sstar_init': init,
                'rounds': rounds, 'stars': len(stars), 'time': time.time() - t1})
    cert_d = {'n': g.n, 'ptr': ptr, 'u': u, 'v': v, 'val': vals, 'b': b,
              'y': np.ones(len(b)), 'bound': float(value)}
    return write_check(cert_d, name, g, tag or 'lsstar', seed, out)


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
        out = run(name, mode, T, seed, tag)
        out.update({'tag': tag, 'commit': commit, 'cpu': cpu_model(),
                    'python': platform.python_version()})
        safe = name.replace('/', '_')
        with open(os.path.join(OUT, f'{tag}_{safe}_{seed}.json'), 'w') as fh:
            json.dump(out, fh)
        print(name, out, flush=True)


if __name__ == '__main__':
    main(sys.argv[4:], float(sys.argv[1]), int(sys.argv[2]), sys.argv[3])

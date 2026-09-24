"""CertiFlip with an independently checked lower-bound certificate.

Runs CertiFlip with budget T, writes the certificate, checks it with
experiments/check_certificate.py (exact arithmetic, row validity) and stores
cost, solver bound and checked bound in one JSON file per run.

usage: colab_run_cert.py T SEED TAG GRAPH..."""
import json, os, platform, sys, time
CC = os.environ.get('CC_ROOT', '/content/cc')
sys.path.insert(0, CC); sys.path.insert(0, os.path.join(CC, 'experiments'))
OUT = os.environ.get('RESULTS_DIR', '/content/results')
import numpy as np


def main(graphs, T, seed, tag):
    import datasets as D
    import ccbench as cc
    from ccbench.certiflip import certiflip
    import check_certificate as C
    os.makedirs(OUT, exist_ok=True)
    try:
        commit = open(os.path.join(CC, 'COMMIT')).read().strip()
    except OSError:
        commit = 'unknown'
    for name in graphs:
        g = D.load(name)
        cert = os.path.join(OUT, f'cert_{name.replace("/", "_")}_{seed}.npz')
        t0 = time.time()
        res = certiflip(g, time_limit=T, rng=seed, cert_path=cert)
        t = time.time() - t0
        assert int(cc.cost(g, res.labels)) == res.cost
        out = {'graph': name, 'n': g.n, 'm': g.m, 'T': T, 'seed': seed, 'tag': tag,
               'commit': commit, 'python': platform.python_version(),
               'cost': int(res.cost), 'lb': float(res.lower_bound), 'time': t}
        if os.path.exists(cert):
            t1 = time.time()
            try:
                rep = C.check(g, dict(np.load(cert)))
                out.update({'check': 'ok', 'certified': int(rep['certified']),
                            'lb_exact': rep['lb_exact'], 'rows': rep['rows'],
                            'star_rows': rep['star_rows'], 'check_time': time.time() - t1})
            except ValueError as exc:
                out.update({'check': f'rejected: {exc}'})
        else:
            out['check'] = 'no bound (support too large)'
        safe = name.replace('/', '_')
        with open(os.path.join(OUT, f'{tag}_{safe}_{seed}.json'), 'w') as fh:
            json.dump(out, fh)
        print(name, out, flush=True)


if __name__ == '__main__':
    main(sys.argv[4:], float(sys.argv[1]), int(sys.argv[2]), sys.argv[3])

"""CertiFlip with an independently checked lower-bound certificate.

Runs CertiFlip with budget T, stamps the certificate with the identity of the
raw instance file (experiments/instance_meta.py), checks it with
experiments/check_certificate.py (which parses the raw file itself), keeps the
certificate under RESULTS_DIR/certs/ for archiving, and stores cost, times,
solver bound and checked bound in one JSON file per run.

usage: colab_run_cert.py T SEED TAG GRAPH..."""
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


def main(graphs, T, seed, tag):
    import datasets as D
    import ccbench as cc
    from ccbench.certiflip import certiflip
    import check_certificate as C
    import instance_meta as M
    os.makedirs(os.path.join(OUT, 'certs'), exist_ok=True)
    try:
        commit = open(os.path.join(CC, 'COMMIT')).read().strip()
    except OSError:
        commit = 'unknown'
    for name in graphs:
        g = D.load(name)
        safe = name.replace('/', '_')
        cert = os.path.join(OUT, f'cert_{safe}_{seed}.npz')
        t0 = time.time()
        res = certiflip(g, time_limit=T, rng=seed, cert_path=cert)
        t = time.time() - t0
        assert int(cc.cost(g, res.labels)) == res.cost
        out = {'graph': name, 'n': g.n, 'm': g.m, 'T': T, 'seed': seed, 'tag': tag,
               'commit': commit, 'python': platform.python_version(), 'cpu': cpu_model(),
               'cost': int(res.cost), 'lb': float(res.lower_bound), 'time': t,
               'lb_time': float(res.lb_time),
               'lp_seed_used': any(h[0] == 'lp-seed' for h in res.history),
               'history': [list(h) for h in res.history]}
        if os.path.exists(cert):
            M.stamp(cert, name, g)
            t1 = time.time()
            try:
                raw = M.raw_file(name)[1]
                rep = C.check_file(raw, cert)
                out.update({'check': 'ok', 'certified': int(rep['certified']),
                            'lb_exact': rep['lb_exact'], 'rows': rep['rows'],
                            'star_rows': rep['star_rows'], 'identity': rep['identity'],
                            'raw_sha256': rep['raw_sha256'], 'edge_sha256': rep['edge_sha256'],
                            'check_time': time.time() - t1})
            except ValueError as exc:
                out.update({'check': f'rejected: {exc}'})
            # archived certificate: downloaded by the fleet manager
            shutil.move(cert, os.path.join(OUT, 'certs', f'{tag}_{safe}_{seed}.npz'))
        else:
            out['check'] = 'no bound (support too large)'
        # the clustering itself, so that the upper bound can be re-checked
        np.savez_compressed(os.path.join(OUT, 'certs', f'{tag}_{safe}_{seed}.labels.npz'),
                            labels=np.asarray(res.labels, dtype=np.int32))
        with open(os.path.join(OUT, f'{tag}_{safe}_{seed}.json'), 'w') as fh:
            json.dump(out, fh)
        print(name, {k: out[k] for k in ('cost', 'lb', 'time', 'check')}, flush=True)


if __name__ == '__main__':
    main(sys.argv[4:], float(sys.argv[1]), int(sys.argv[2]), sys.argv[3])

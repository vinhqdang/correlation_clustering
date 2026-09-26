"""Head-to-head on one machine, one solver at a time.

KaPoCE and PXMem run one after the other, each pinned to CPU 0 while the rest
of the machine is idle, with the same wall-clock budget.  KaPoCE is stopped as
in PACE 2021: SIGTERM at the limit, after which its output is read to the end
by a reader that is never starved.  Its seed is taken from the run seed (the
only change to its source, see ``prepare``).  The output is checked: the edit
set is valid iff the edited graph is a disjoint union of cliques, i.e. iff the
number of edits equals the cost of the clustering it induces.

Protocol v2 (tags starting with "s2"; see Section 5 of the paper):
  * both solvers read the same PACE-format text file inside their measured time
    (KaPoCE from stdin, PXMem parses it with ccbench.read_pace);
  * KaPoCE's internal time limits are set from T (590/600 and 90/600 of T, the
    ratios of its PACE submission) and it gets SIGTERM at T;
  * the order of the two solvers is random (fixed by graph and seed);
  * PXMem is scored at KaPoCE's elapsed wall-clock time, from its full
    improvement history, and also at T;
  * wall-clock and CPU time (wait4) are recorded for both.

usage: colab_run_seq.py T SEED TAG GRAPH...
options can be appended to the tag as "+key=value" and are passed to pxmem."""
import json, os, platform, subprocess, sys, tempfile, threading, time
CC = os.environ.get('CC_ROOT', '/content/cc')
sys.path.insert(0, CC); sys.path.insert(0, os.path.join(CC, 'experiments'))

KSRC = os.environ.get('KAPOCE_SRC', '/content/kapoce')
KBIN = os.path.join(KSRC, 'build_seed', 'ClusterEditing')
KBIN2 = os.path.join(KSRC, 'build_v2', 'ClusterEditing')
PACE_DIR = os.path.join(CC, 'data', 'raw', 'pace')
OUT = os.environ.get('RESULTS_DIR', '/content/results')
READY = os.environ.get('SEQ_READY', '/content/SEQ_READY2')
CPU = '0'
GRACE = 60  # seconds KaPoCE may take to write its solution after SIGTERM


def sh(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True)


def prepare():
    """Seeded KaPoCE build and PACE heuristic instances, once per machine."""
    lock = READY
    if os.path.exists(lock):
        return
    src = os.path.join(KSRC, 'cluster_editing', 'application', 'heuristic.cc')
    if 'KAPOCE_SEED' not in open(src).read():
        sh(f"sed -i -e 's/context.general.seed = 1;/context.general.seed = "
           f"getenv(\"KAPOCE_SEED\") ? atoi(getenv(\"KAPOCE_SEED\")) : 1;/' "
           f"-e '1i #include <cstdlib>' {src}")
    if not os.path.exists(KBIN):
        r = sh(f"mkdir -p {KSRC}/build_seed && cd {KSRC}/build_seed && cmake .. "
               "-DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_CXX_FLAGS=\"-include cstdint -include limits "
               "-include string -include memory -include algorithm -include functional "
               "-include stdexcept\" > cmake.log 2>&1 && make -j2 ClusterEditing > make.log 2>&1")
        if not os.path.exists(KBIN):
            raise RuntimeError('seeded KaPoCE build failed: ' + r.stderr[-500:])
    if not os.path.exists(KBIN2):
        # v2: time limits from the environment (defaults: the PACE submission)
        sh(f"sed -i -e 's/context.general.time_limit = 590;/context.general.time_limit = "
           f"getenv(\"KAPOCE_TIME\") ? atof(getenv(\"KAPOCE_TIME\")) : 590;/' "
           f"-e 's/context.refinement.evo.time_limit = 90;/context.refinement.evo.time_limit = "
           f"getenv(\"KAPOCE_EVO\") ? atof(getenv(\"KAPOCE_EVO\")) : 90;/' {src}")
        r = sh(f"mkdir -p {KSRC}/build_v2 && cd {KSRC}/build_v2 && cmake .. "
               "-DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_CXX_FLAGS=\"-include cstdint -include limits "
               "-include string -include memory -include algorithm -include functional "
               "-include stdexcept\" > cmake.log 2>&1 && make -j2 ClusterEditing > make.log 2>&1")
        if not os.path.exists(KBIN2):
            raise RuntimeError('KaPoCE v2 build failed: ' + r.stderr[-500:])
    if not os.path.isdir(os.path.join(PACE_DIR, 'heur')):
        os.makedirs(PACE_DIR, exist_ok=True)
        sh('cd /tmp && rm -rf pace21 && git clone -q --depth 1 '
           'https://github.com/PACE-challenge/Cluster-Editing-PACE-2021-instances pace21')
        for track in ('exact', 'heur'):
            d = os.path.join(PACE_DIR, track)
            os.makedirs(d, exist_ok=True)
            sh(f"cd /tmp/pace21 && for f in $(find . -name '{track}*.gr*'); do cp $f {d}/; done; "
               f"cd {d} && for f in *.gz; do [ -e \"$f\" ] && gunzip -f \"$f\"; done; "
               f"for f in *.xz; do [ -e \"$f\" ] && unxz -f \"$f\"; done")
    # compile the Numba kernels once, outside any timed run
    sh(f"cd {CC} && taskset -c {CPU} python3 -c \"import ccbench as cc; "
       "from ccbench.memetic import pxmem; from ccbench.generators import planted_partition; "
       "g = planted_partition(300, 10, 0.6, 0.02, rng=0)[0]; pxmem(g, 5, rng=0)\"")
    open(lock, 'w').close()


def run_kapoce(g, T, seed):
    import baselines as B
    with tempfile.TemporaryDirectory() as d:
        inp = os.path.join(d, 'in.gr')
        B.write_pace(g, inp)
        env = dict(os.environ, KAPOCE_SEED=str(seed + 1))
        t0 = time.time()
        with open(inp) as fin:
            p = subprocess.Popen(['taskset', '-c', CPU, KBIN], stdin=fin, stdout=subprocess.PIPE,
                                 stderr=subprocess.DEVNULL, env=env)
        chunks = []
        reader = threading.Thread(target=lambda: chunks.append(p.stdout.read()))
        reader.start()
        stopped = False

        def reap(deadline):
            # wait4 instead of Popen.wait, to get the resource usage of the child
            while time.time() < deadline:
                pid, status, ru = os.wait4(p.pid, os.WNOHANG)
                if pid:
                    return status, ru
                time.sleep(0.05)
            return None

        done = reap(t0 + T)
        if done is None:
            p.terminate()
            stopped = True
            done = reap(time.time() + GRACE)
            if done is None:
                p.kill()
                _, status, ru = os.wait4(p.pid, 0)
                done = status, ru
        status, ru = done
        p.returncode = status  # keep Popen from waiting for the reaped child
        reader.join()
        elapsed = time.time() - t0
    rows = []
    for line in chunks[0].decode().splitlines():
        parts = line.split()
        if len(parts) == 2 and parts[0].isdigit() and parts[1].isdigit():
            rows.append((int(parts[0]) - 1, int(parts[1]) - 1))
    import ccbench as cc
    lab = B.apply_edits(g, rows)
    c = int(cc.cost(g, lab))
    n_edits = len({(min(u, v), max(u, v)) for u, v in rows})
    return {'kapoce': c, 'kapoce_time': elapsed, 'kapoce_sigterm': stopped,
            'kapoce_valid': bool(n_edits == c),
            'kapoce_edits': n_edits, 'kapoce_rss_mb': ru.ru_maxrss / 1024,
            'kapoce_exit': status}


PX_CHILD = r'''
import json, resource, sys, time
sys.path.insert(0, a_cc := sys.argv[2]); sys.path.insert(0, a_cc + '/experiments')
import datasets as D, ccbench as cc
from ccbench.memetic import pxmem
a = json.loads(sys.argv[1])
g = D.load(a['graph'])
t0 = time.time()
hist = []
lab, c = pxmem(g, a['T'], rng=a['seed'], history=hist, **a['opts'])
t = time.time() - t0
assert int(cc.cost(g, lab)) == c
print('RESULT' + json.dumps({'ours': int(c), 'ours_time': t,
      'ours_rss_mb': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024,
      'hist': hist[::max(1, len(hist) // 50)]}))
'''


def run_ours(name, T, seed, opts):
    arg = json.dumps({'graph': name, 'T': T, 'seed': seed, 'opts': opts})
    r = subprocess.run(['taskset', '-c', CPU, sys.executable, '-c', PX_CHILD, arg, CC],
                       capture_output=True, text=True, cwd=CC)
    for line in r.stdout.splitlines():
        if line.startswith('RESULT'):
            return json.loads(line[6:])
    raise RuntimeError('pxmem failed: ' + r.stderr[-2000:])


def machine():
    lscpu = sh('lscpu').stdout
    keep = ('Model name', 'Thread(s) per core', 'Core(s) per socket', 'Socket(s)', 'CPU MHz',
            'L3 cache')
    info = {l.split(':')[0].strip(): l.split(':', 1)[1].strip() for l in lscpu.splitlines()
            if ':' in l and l.split(':')[0].strip() in keep}
    info['mem_gb'] = round(os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES') / 2**30, 1)
    info['python'] = platform.python_version()
    return info


def main(graphs, T, seed, tag):
    import datasets as D
    prepare()
    os.makedirs(OUT, exist_ok=True)
    opts = {}
    for part in tag.split('+')[1:]:
        k, v = part.split('=')
        v = float(v)
        opts[k] = int(v) if v.is_integer() else v
    try:
        commit = open(os.path.join(CC, 'COMMIT')).read().strip()
    except OSError:
        commit = 'unknown'
    for name in graphs:
        g = D.load(name)
        res = {'graph': name, 'n': g.n, 'm': g.m, 'T': T, 'seed': seed, 'tag': tag,
               'opts': opts, 'protocol': 'sequential, pinned to one CPU', 'commit': commit,
               'machine': machine()}
        res.update(run_kapoce(g, T, seed))
        res.update(run_ours(name, T, seed, opts))
        safe = name.replace('/', '_')
        with open(os.path.join(OUT, f'{tag}_{safe}_{seed}.json'), 'w') as fh:
            json.dump(res, fh)
        print(name, 'ours', res['ours'], 'kapoce', res['kapoce'], res['kapoce_valid'], flush=True)


# ---------------------------------------------------------------------------
# protocol v2

def _wait(p, deadline):
    while time.time() < deadline:
        pid, status, ru = os.wait4(p.pid, os.WNOHANG)
        if pid:
            return status, ru
        time.sleep(0.02)
    return None


def run_kapoce2(g, inp, T, seed):
    """KaPoCE reading the .gr file from stdin; time limits set from T."""
    import baselines as B
    import ccbench as cc
    env = dict(os.environ, KAPOCE_SEED=str(seed + 1), KAPOCE_TIME=repr(T * 590 / 600),
               KAPOCE_EVO=repr(T * 90 / 600))
    t0 = time.time()
    with open(inp) as fin:
        p = subprocess.Popen(['taskset', '-c', CPU, KBIN2], stdin=fin, stdout=subprocess.PIPE,
                             stderr=subprocess.DEVNULL, env=env)
    chunks = []
    reader = threading.Thread(target=lambda: chunks.append(p.stdout.read()))
    reader.start()
    stopped = False
    done = _wait(p, t0 + T)
    if done is None:
        p.terminate()
        stopped = True
        done = _wait(p, time.time() + GRACE)
        if done is None:
            p.kill()
            _, status, ru = os.wait4(p.pid, 0)
            done = status, ru
    status, ru = done
    p.returncode = status
    reader.join()
    elapsed = time.time() - t0
    rows = []
    for line in chunks[0].decode().splitlines():
        parts = line.split()
        if len(parts) == 2 and parts[0].isdigit() and parts[1].isdigit():
            rows.append((int(parts[0]) - 1, int(parts[1]) - 1))
    lab = B.apply_edits(g, rows)
    c = int(cc.cost(g, lab))
    n_edits = len({(min(u, v), max(u, v)) for u, v in rows})
    return {'_labels_kapoce': lab,
            'kapoce': c, 'kapoce_time': elapsed, 'kapoce_cpu': ru.ru_utime + ru.ru_stime,
            'kapoce_sigterm': stopped, 'kapoce_valid': bool(n_edits == c),
            'kapoce_edits': n_edits, 'kapoce_rss_mb': ru.ru_maxrss / 1024,
            'kapoce_exit': status}


PX_CHILD2 = r"""
import json, sys, time
t_start = float(sys.argv[3])
sys.path.insert(0, a_cc := sys.argv[2]); sys.path.insert(0, a_cc + '/experiments')
import ccbench as cc
from ccbench.memetic import pxmem
a = json.loads(sys.argv[1])
g = cc.read_pace(a['input'])
off = time.time() - t_start
hist = []
lab, c = pxmem(g, max(1.0, a['T'] - off), rng=a['seed'], history=hist, **a['opts'])
t = time.time() - t_start
assert int(cc.cost(g, lab)) == c
if a.get('save'):
    import numpy as np
    np.savez_compressed(a['save'], labels=np.asarray(lab, dtype=np.int32))
print('RESULT' + json.dumps({'ours_final': int(c), 'ours_time': t, 'ours_load': off,
      'hist': [(off + h[0], int(h[1])) for h in hist]}))
"""


def run_ours2(inp, T, seed, opts, save=None):
    arg = json.dumps({'input': inp, 'T': T, 'seed': seed, 'opts': opts, 'save': save})
    t0 = time.time()
    p = subprocess.Popen(['taskset', '-c', CPU, sys.executable, '-c', PX_CHILD2, arg, CC,
                          repr(t0)], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         text=True, cwd=CC)
    out, err = [], []
    tr = [threading.Thread(target=lambda: out.append(p.stdout.read())),
          threading.Thread(target=lambda: err.append(p.stderr.read()))]
    for x in tr:
        x.start()
    # PXMem stops itself at T; the hard limit only guards against a hang
    done = _wait(p, t0 + T + 600)
    if done is None:
        p.kill()
        done = os.wait4(p.pid, 0)[1:]
    status, ru = done
    p.returncode = status
    for x in tr:
        x.join()
    for line in out[0].splitlines():
        if line.startswith('RESULT'):
            r = json.loads(line[6:])
            r['ours_cpu'] = ru.ru_utime + ru.ru_stime
            r['ours_rss_mb'] = ru.ru_maxrss / 1024
            return r
    raise RuntimeError('pxmem failed: ' + err[0][-2000:])


def at(hist, t):
    """Best cost found by time t (None if there was no solution yet)."""
    best = None
    for tt, c in hist:
        if tt <= t:
            best = c if best is None else min(best, c)
    return best


def main2(graphs, T, seed, tag, opts):
    import datasets as D
    import baselines as B
    import hashlib
    prepare()
    os.makedirs(OUT, exist_ok=True)
    try:
        commit = open(os.path.join(CC, 'COMMIT')).read().strip()
    except OSError:
        commit = 'unknown'
    for name in graphs:
        g = D.load(name)
        safe = name.replace('/', '_')
        inp = os.path.join('/tmp', f'in_{safe}.gr')
        B.write_pace(g, inp)
        kfirst = int(hashlib.sha256(f'{name}:{seed}'.encode()).hexdigest(), 16) % 2 == 0
        res = {'graph': name, 'n': g.n, 'm': g.m, 'T': T, 'seed': seed, 'tag': tag,
               'opts': opts, 'protocol': 'v2: sequential, pinned to one CPU, same text input, '
               'KaPoCE limits from T, random order, PXMem scored at KaPoCE time',
               'kapoce_first': kfirst, 'commit': commit, 'machine': machine()}
        # the final clusterings of the primary runs (seed 0 at 600 s on the SNAP
        # graphs) are archived so that their costs can be re-checked
        keep = seed == 0 and tag == 's2h'
        save = f'/tmp/px_{safe}.npz' if keep else None
        if kfirst:
            res.update(run_kapoce2(g, inp, T, seed))
            res.update(run_ours2(inp, T, seed, opts, save))
        else:
            res.update(run_ours2(inp, T, seed, opts, save))
            res.update(run_kapoce2(g, inp, T, seed))
        lab_k = res.pop('_labels_kapoce')
        if keep:
            import instance_meta as M
            import numpy as np
            meta = {k: np.array(v) for k, v in M.meta(name, g).items()}
            os.makedirs(os.path.join(OUT, 'certs'), exist_ok=True)
            np.savez_compressed(os.path.join(OUT, 'certs', f'{tag}_{safe}_{seed}.kapoce.labels.npz'),
                                labels=np.asarray(lab_k, dtype=np.int32), **meta)
            np.savez_compressed(os.path.join(OUT, 'certs', f'{tag}_{safe}_{seed}.pxmem.labels.npz'),
                                labels=np.load(save)['labels'], **meta)
            os.remove(save)
        cut = min(T, res['kapoce_time'])
        res['cut'] = cut
        res['ours'] = at(res['hist'], cut)          # the reported comparison
        res['ours_at_T'] = at(res['hist'], T)
        os.remove(inp)
        with open(os.path.join(OUT, f'{tag}_{safe}_{seed}.json'), 'w') as fh:
            json.dump(res, fh)
        print(name, 'ours', res['ours'], 'kapoce', res['kapoce'], res['kapoce_valid'], flush=True)


if __name__ == '__main__':
    T_, seed_, tag_, graphs_ = float(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4:]
    if tag_.startswith('s2'):
        opts_ = {}
        for part in tag_.split('+')[1:]:
            k_, v_ = part.split('=')
            v_ = float(v_)
            opts_[k_] = int(v_) if v_.is_integer() else v_
        main2(graphs_, T_, seed_, tag_, opts_)
    else:
        main(graphs_, T_, seed_, tag_)

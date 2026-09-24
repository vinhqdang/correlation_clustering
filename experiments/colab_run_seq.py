"""Head-to-head on one machine, one solver at a time.

KaPoCE and PXMem run one after the other, each pinned to CPU 0 while the rest
of the machine is idle, with the same wall-clock budget.  KaPoCE is stopped as
in PACE 2021: SIGTERM at the limit, after which its output is read to the end
by a reader that is never starved.  Its seed is taken from the run seed (the
only change to its source, see ``prepare``).  The output is checked: the edit
set is valid iff the edited graph is a disjoint union of cliques, i.e. iff the
number of edits equals the cost of the clustering it induces.

usage: colab_run_seq.py T SEED TAG GRAPH...
options can be appended to the tag as "+key=value" and are passed to pxmem."""
import json, os, platform, subprocess, sys, tempfile, threading, time
CC = os.environ.get('CC_ROOT', '/content/cc')
sys.path.insert(0, CC); sys.path.insert(0, os.path.join(CC, 'experiments'))

KSRC = os.environ.get('KAPOCE_SRC', '/content/kapoce')
KBIN = os.path.join(KSRC, 'build_seed', 'ClusterEditing')
PACE_DIR = os.path.join(CC, 'data', 'raw', 'pace')
OUT = os.environ.get('RESULTS_DIR', '/content/results')
READY = os.environ.get('SEQ_READY', '/content/SEQ_READY')
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


if __name__ == '__main__':
    main(sys.argv[4:], float(sys.argv[1]), int(sys.argv[2]), sys.argv[3])

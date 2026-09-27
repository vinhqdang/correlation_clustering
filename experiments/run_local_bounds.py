"""Bound jobs of the PACE exact track on a local multi-core machine, several
single-threaded runs in parallel (the same runner as on the Colab machines).

  python experiments/run_local_bounds.py [--jobs N] [--only lwarm|ltri]

Jobs: on every fourth instance (exact004, exact008, ..., exact200)
  lwarm  the CertiFlip bound with 1800 s and one block covering the graph
         (checked certificate);
  ltri   the metric LP on P with triangle rows, 900 s (LP value).
Results go to results/local/ (JSON, certificates under certs/); a job whose
result exists there or in results/colab/ is skipped, so the script can be
restarted.  Commit and push results/local afterwards.
"""
import argparse
import os
import platform
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, "results", "local")
JOBS = [("lwarm", 1800), ("ltri", 900)]


def result_path(d, tag, graph, seed=0):
    return os.path.join(d, f"{tag}_{graph.replace('/', '_')}_{seed}.json")


def run(tag, T, graph):
    env = dict(os.environ, CC_ROOT=ROOT, RESULTS_DIR=OUT, OMP_NUM_THREADS="1",
               NUMBA_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    log = os.path.join(OUT, "logs", f"{tag}_{graph.replace('/', '_')}.txt")
    t0 = time.time()
    with open(log, "w") as fh:
        p = subprocess.run([sys.executable, os.path.join(ROOT, "experiments", "colab_run_lp.py"),
                            str(T), "0", tag, graph], stdout=fh, stderr=subprocess.STDOUT, env=env)
    return tag, graph, p.returncode, time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) // 2),
                    help="parallel runs (default: half of the logical cores)")
    ap.add_argument("--only", choices=[t for t, _ in JOBS])
    ap.add_argument("--limit", type=int, default=0, help="run at most this many jobs (a test)")
    a = ap.parse_args()
    os.makedirs(os.path.join(OUT, "logs"), exist_ok=True)
    graphs = [f"pace-exact/exact{i:03d}.gr" for i in range(4, 201, 4)]
    todo = []
    for tag, T in JOBS:
        if a.only and tag != a.only:
            continue
        for g in graphs:
            if any(os.path.exists(result_path(d, tag, g))
                   for d in (OUT, os.path.join(ROOT, "results", "colab"))):
                continue
            todo.append((tag, T, g))
    if a.limit:
        todo = todo[:a.limit]
    # the PACE instances are fetched once, before the parallel runs
    sys.path[:0] = [ROOT, os.path.join(ROOT, "experiments")]
    os.environ["CC_ROOT"] = ROOT
    import colab_run_lp
    colab_run_lp.ensure_pace()
    est = sum(T for _, T, _ in todo) / max(1, a.jobs) / 3600
    print(f"{len(todo)} jobs, {a.jobs} in parallel, about {est:.1f} h "
          f"({platform.node()}, {os.cpu_count()} logical cores)", flush=True)
    done = 0
    with ProcessPoolExecutor(max_workers=a.jobs) as ex:
        futs = [ex.submit(run, *j) for j in todo]
        for f in as_completed(futs):
            tag, g, rc, dt = f.result()
            done += 1
            ok = os.path.exists(result_path(OUT, tag, g))
            print(f"[{time.strftime('%H:%M:%S')}] {done}/{len(todo)} {tag} {g}: "
                  f"{'ok' if ok and rc == 0 else f'FAILED (exit {rc}), see results/local/logs'}"
                  f" ({dt:.0f} s)", flush=True)
    print("done; now: git add results/local && git commit -m 'PACE bound runs' && git push")


if __name__ == "__main__":
    main()

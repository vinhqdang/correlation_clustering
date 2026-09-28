"""Star local search with the stopping rule of the KaPoCE root bound (tag
lsstar+kr) on the PACE exact track, seeds 0-4, on this machine; skips runs
whose result file exists (results/colab or results/local).  The rule counts
rounds, so the result does not depend on the machine unless the 60 s cap is
reached."""
import os
import subprocess
import sys
from multiprocessing import Pool

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
TAG, T = "lsstar+kr", 60


def one(job):
    seed, inst = job
    name = f"pace-exact/{inst}"
    stem = f"{TAG}_{name.replace('/', '_')}_{seed}.json"
    if any(os.path.exists(os.path.join(ROOT, "results", d, stem)) for d in ("colab", "local")):
        return
    env = dict(os.environ, CC_ROOT=ROOT, RESULTS_DIR=os.path.join(ROOT, "results", "local"))
    subprocess.run([sys.executable, os.path.join(ROOT, "experiments", "colab_run_lp.py"), str(T),
                    str(seed), TAG, name], env=env, cwd=ROOT, capture_output=True, timeout=3600)


if __name__ == "__main__":
    insts = sorted(os.listdir(os.path.join(ROOT, "data", "raw", "_pace21", "exact")))
    jobs = [(s, i) for s in range(5) for i in insts]
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 3) as p:
        for _ in p.imap_unordered(one, jobs):
            pass

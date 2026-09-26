"""SCC (Hausberger, Faraj, Schulz, ALENEX 2025) on complete instances.

SCC optimises correlation clustering on sparse signed graphs, where a pair
without an edge costs nothing.  To make it optimise the complete-graph
objective, every non-adjacent pair is written as an explicit edge of weight -1
(positive edges have weight +1).  This needs n(n-1)/2 edges and is therefore
done only on small instances.  The cost of SCC's clustering is recomputed by
the independent evaluator of check_certificate.py.

usage: run_scc.py T SEED MODE GRAPH...     (MODE: evo = scc_evolutionary with
       one MPI process and time limit T;  ml = one multilevel run of scc)"""
import json
import os
import subprocess
import sys
import time

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path[:0] = [ROOT, os.path.join(ROOT, "experiments")]
SCC = os.environ.get("SCC_DIR", "/tmp/scc/build")
OUT = os.path.join(ROOT, "results", "scc")


def write_complete(g, path):
    n = g.n
    with open(path, "w") as fh:
        fh.write(f"{n} {n * (n - 1) // 2} 1\n")
        for u in range(n):
            w = -np.ones(n, dtype=np.int64)
            w[g.indices[g.indptr[u]:g.indptr[u + 1]]] = 1
            w[u] = 0
            idx = np.flatnonzero(w)
            fh.write(" ".join(f"{v + 1} {w[v]}" for v in idx.tolist()) + "\n")


def main(T, seed, mode, graphs):
    import datasets as D
    import instance_meta as M
    import check_certificate as C
    os.makedirs(OUT, exist_ok=True)
    for name in graphs:
        safe = name.replace("/", "_")
        out = os.path.join(OUT, f"scc{mode}_{safe}_{seed}.json")
        if os.path.exists(out):
            continue
        g = D.load(name)
        gpath = f"/tmp/scc_{safe}.graph"
        cpath = f"/tmp/scc_{safe}_{seed}.clu"
        write_complete(g, gpath)
        if mode == "evo":
            cmd = ["mpirun", "--allow-run-as-root", "--oversubscribe", "-n", "1",
                   os.path.join(SCC, "scc_evolutionary_int"), gpath, f"--seed={seed}",
                   f"--time_limit={T}", f"--output_filename={cpath}"]
        else:
            cmd = [os.path.join(SCC, "scc_int"), gpath, f"--seed={seed}",
                   f"--output_filename={cpath}"]
        t0 = time.time()
        r = subprocess.run(cmd, capture_output=True, text=True)
        t = time.time() - t0
        os.remove(gpath)
        rec = {"graph": name, "n": g.n, "m": g.m, "mode": mode, "T": T, "seed": seed,
               "time": t, "exit": r.returncode}
        for line in r.stdout.splitlines():
            if line.startswith("cut"):
                # SCC's objective: weight of the cut pairs; cost = cut + #negative pairs
                rec["scc_cost"] = int(float(line.split()[1])) + g.n * (g.n - 1) // 2 - g.m
        if r.returncode == 0 and os.path.exists(cpath):
            lab = np.loadtxt(cpath, dtype=np.int64).reshape(-1)
            rel, raw, fmt, sc = M.raw_file(name)
            rec["cost"] = int(C.clustering_cost(C.read_instance(raw, fmt, sc), lab))
            os.remove(cpath)
        else:
            rec["stderr"] = r.stderr[-2000:]
        json.dump(rec, open(out, "w"))
        print(name, {k: rec.get(k) for k in ("cost", "time", "exit")}, flush=True)


if __name__ == "__main__":
    main(float(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4:])

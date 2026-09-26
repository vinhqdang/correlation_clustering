"""Dataset registry for the benchmark."""
from __future__ import annotations

import os
import numpy as np

import ccbench as cc
from ccbench.generators import planted_partition, sparse_planted

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "raw")
SNAP = "https://snap.stanford.edu/data/"

REAL = {
    # name: (file, url suffix, reader kwargs)
    "ca-GrQc": ("ca-GrQc.txt.gz", "ca-GrQc.txt.gz", {}),
    "ca-HepTh": ("ca-HepTh.txt.gz", "ca-HepTh.txt.gz", {}),
    "ca-HepPh": ("ca-HepPh.txt.gz", "ca-HepPh.txt.gz", {}),
    "ca-AstroPh": ("ca-AstroPh.txt.gz", "ca-AstroPh.txt.gz", {}),
    "ca-CondMat": ("ca-CondMat.txt.gz", "ca-CondMat.txt.gz", {}),
    "email-Enron": ("email-Enron.txt.gz", "email-Enron.txt.gz", {}),
    "loc-Brightkite": ("loc-brightkite_edges.txt.gz", "loc-brightkite_edges.txt.gz", {}),
    "soc-Epinions": ("soc-Epinions1.txt.gz", "soc-Epinions1.txt.gz", {}),
    "Slashdot+": ("soc-sign-Slashdot090221.txt.gz", "soc-sign-Slashdot090221.txt.gz",
                  {"sign_col": 2}),
    "Epinions+": ("soc-sign-epinions.txt.gz", "soc-sign-epinions.txt.gz", {"sign_col": 2}),
    "BitcoinOTC+": ("soc-sign-bitcoinotc.csv.gz", "soc-sign-bitcoinotc.csv.gz", {"sign_col": 2}),
    "BitcoinAlpha+": ("soc-sign-bitcoinalpha.csv.gz", "soc-sign-bitcoinalpha.csv.gz",
                      {"sign_col": 2}),
    "com-Amazon": ("com-amazon.ungraph.txt.gz", "bigdata/communities/com-amazon.ungraph.txt.gz", {}),
    "com-DBLP": ("com-dblp.ungraph.txt.gz", "bigdata/communities/com-dblp.ungraph.txt.gz", {}),
    "com-Youtube": ("com-youtube.ungraph.txt.gz", "bigdata/communities/com-youtube.ungraph.txt.gz",
                    {}),
    # held-out graphs for the head-to-head (fixed before any run, see experiments/heldout.md)
    "email-Eu-core": ("email-Eu-core.txt.gz", "email-Eu-core.txt.gz", {}),
    "facebook": ("facebook_combined.txt.gz", "facebook_combined.txt.gz", {}),
    "wiki-Vote": ("wiki-Vote.txt.gz", "wiki-Vote.txt.gz", {}),
    "cit-HepTh": ("cit-HepTh.txt.gz", "cit-HepTh.txt.gz", {}),
    "cit-HepPh": ("cit-HepPh.txt.gz", "cit-HepPh.txt.gz", {}),
    "p2p-Gnutella31": ("p2p-Gnutella31.txt.gz", "p2p-Gnutella31.txt.gz", {}),
    "soc-Slashdot0902": ("soc-Slashdot0902.txt.gz", "soc-Slashdot0902.txt.gz", {}),
    "loc-Gowalla": ("loc-gowalla_edges.txt.gz", "loc-gowalla_edges.txt.gz", {}),
    "email-EuAll": ("email-EuAll.txt.gz", "email-EuAll.txt.gz", {}),
    "amazon0302": ("amazon0302.txt.gz", "amazon0302.txt.gz", {}),
    "web-NotreDame": ("web-NotreDame.txt.gz", "web-NotreDame.txt.gz", {}),
    "roadNet-PA": ("roadNet-PA.txt.gz", "roadNet-PA.txt.gz", {}),
}


def ensure(name: str) -> str:
    f, url, _ = REAL[name]
    path = os.path.join(ROOT, f)
    if not os.path.exists(path):
        import urllib.request
        os.makedirs(ROOT, exist_ok=True)
        urllib.request.urlretrieve(SNAP + url, path)
    return path


def load(name: str) -> cc.Graph:
    if name in REAL:
        path = ensure(name)
        cache = path + ".npz"
        if os.path.exists(cache):
            z = np.load(cache)
            return cc.Graph(int(z["n"]), z["indptr"], z["indices"])
        g = cc.read_edgelist(path, **REAL[name][2])
        np.savez(cache, n=g.n, indptr=g.indptr, indices=g.indices)
        return g
    if name.startswith("pace-"):
        return cc.read_pace(os.path.join(ROOT, "pace", name[5:]))
    return synthetic(name)[0]


def synthetic(name: str):
    """Names: sbm-<n>-<k>-<pin>-<pout>-<seed>, sparse-<n>-<avg>-<degout>-<pin>-<seed>."""
    parts = name.split("-")
    kind = parts[0]
    if kind == "sbm":
        n, k = int(parts[1]), int(parts[2])
        pin, pout, seed = float(parts[3]), float(parts[4]), int(parts[5])
        return planted_partition(0, 0, pin, pout, rng=seed, sizes=np.full(k, n // k))
    if kind == "sparse":
        n, avg = int(parts[1]), int(parts[2])
        dout, pin, seed = float(parts[3]), float(parts[4]), int(parts[5])
        return sparse_planted(n, avg, dout, pin, rng=seed)
    raise KeyError(name)

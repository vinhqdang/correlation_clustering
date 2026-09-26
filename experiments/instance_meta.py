"""Identity of an instance, recorded in every certificate.

The certificate stores the raw file it was computed for (name, SHA-256,
format) together with n, m and a hash of the solver's own edge list.  The
checker parses the raw file itself and must obtain the same values."""
from __future__ import annotations

import hashlib
import os

import numpy as np

import datasets as D


def edge_hash(n, edges):
    """SHA-256 of n and the lexicographically sorted edge list (u < v, int64)."""
    e = np.asarray(edges, dtype=np.int64).reshape(-1, 2)
    e = np.sort(e, axis=1)
    e = e[np.lexsort((e[:, 1], e[:, 0]))]
    h = hashlib.sha256()
    h.update(np.int64(n).tobytes())
    h.update(np.ascontiguousarray(e, dtype="<i8").tobytes())
    return h.hexdigest()


def raw_file(name: str):
    """(path relative to data/raw, absolute path, format, sign column)."""
    if name in D.REAL:
        f, _, kw = D.REAL[name]
        return f, D.ensure(name), "snap", int(kw.get("sign_col", -1))
    if name.startswith("pace-"):
        track, f = name[5:].split("/")
        rel = os.path.join("pace", track, f)
        return rel, os.path.join(D.ROOT, rel), "pace", -1
    raise KeyError(f"{name} has no raw file")


def meta(name: str, g) -> dict:
    rel, path, fmt, sign_col = raw_file(name)
    with open(path, "rb") as fh:
        sha = hashlib.sha256(fh.read()).hexdigest()
    return {"instance_name": name, "instance_file": rel, "instance_format": fmt,
            "sign_col": sign_col, "raw_sha256": sha, "instance_m": int(g.m),
            "edge_sha256": edge_hash(g.n, g.edges())}


def stamp(cert_path: str, name: str, g) -> None:
    """Add the instance identity to a certificate written by CertiFlip."""
    cert = dict(np.load(cert_path))
    cert.update({k: np.array(v) for k, v in meta(name, g).items()})
    np.savez_compressed(cert_path, **cert)

"""Download every instance used in the paper and check it against
data/MANIFEST.sha256 (SHA-256 of the raw files as distributed).

SNAP graphs come from https://snap.stanford.edu/data/, the PACE 2021 instances
from the challenge repository at a fixed commit.

usage: python experiments/fetch_data.py [--write-manifest]"""
import hashlib
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import datasets as D  # noqa: E402

PACE_REPO = "https://github.com/PACE-challenge/Cluster-Editing-PACE-2021-instances"
PACE_COMMIT = "67f5df368c400f15fa2527ea1eb48ac1bd6a2d7a"
MANIFEST = os.path.join(HERE, "..", "data", "MANIFEST.sha256")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch_pace():
    d = os.path.join(D.ROOT, "pace")
    if all(os.path.isdir(os.path.join(d, t)) and len(os.listdir(os.path.join(d, t))) >= 200
           for t in ("exact", "heur")):
        return
    tmp = os.path.join(D.ROOT, "_pace21")
    if not os.path.isdir(tmp):
        subprocess.run(["git", "clone", "-q", PACE_REPO, tmp], check=True)
    subprocess.run(["git", "-C", tmp, "checkout", "-q", PACE_COMMIT], check=True)
    for track in ("exact", "heur"):
        t = os.path.join(d, track)
        os.makedirs(t, exist_ok=True)
        subprocess.run(f"for f in $(find {tmp} -name '{track}*.gr*'); do cp $f {t}/; done; "
                       f"cd {t} && for f in *.gz; do [ ! -e \"$f\" ] || gunzip -f \"$f\"; done; "
                       f"for f in *.xz; do [ ! -e \"$f\" ] || unxz -f \"$f\"; done",
                       shell=True, check=True)


def files():
    out = [D.REAL[name][0] for name in D.REAL]
    for track in ("exact", "heur"):
        out += [os.path.join("pace", track, f"{track}{i:03d}.gr") for i in range(1, 201)]
    return out


def main(write):
    for name in D.REAL:
        D.ensure(name)
    fetch_pace()
    have = {f: sha256(os.path.join(D.ROOT, f)) for f in files()}
    if write:
        with open(MANIFEST, "w") as fh:
            for f in sorted(have):
                fh.write(f"{have[f]}  {f}\n")
        print(f"wrote {len(have)} hashes to {MANIFEST}")
        return 0
    want = dict(line.split()[::-1] for line in open(MANIFEST) if line.strip())
    bad = [f for f in want if have.get(f) != want[f]]
    for f in bad:
        print(f"MISMATCH {f}: {have.get(f)} != {want[f]}")
    print(f"{len(want) - len(bad)} of {len(want)} instance files match the manifest")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main("--write-manifest" in sys.argv))

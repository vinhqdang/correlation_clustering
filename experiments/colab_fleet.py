"""Self-healing job runner on a fleet of Google Colab CPU sessions.

Every POLL seconds, for every session:
  * probe it (setup finished? jobs running? result files present?);
  * download finished result files to results/colab/;
  * launch the next pending job if the machine is idle;
  * if the probe fails repeatedly, recreate the session, re-run the setup and
    put its unfinished jobs back into the queue.
New results are committed and pushed periodically.  The queue lives in
results/colab/queue.json and can be edited while the fleet is running
(add jobs with ``python experiments/colab_fleet.py add ...``).

Run under the supervisor script experiments/colab_fleet.sh, which restarts
this process if it exits.
"""
from __future__ import annotations

import json
import os
from concurrent.futures import ThreadPoolExecutor
import subprocess
import sys
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, "results", "colab")
QUEUE = os.path.join(OUT, "queue.json")
STATE = os.path.join(OUT, "fleet_state.json")
INBOX = os.path.join(OUT, "inbox")  # job batches added while the fleet runs
LOG = os.path.join(OUT, "fleet.log")
COLAB = os.path.expanduser("/root/.local/bin/colab")
COLAB_PY = "/root/.local/share/uv/tools/google-colab-cli/bin/python"
REAP = os.path.join(ROOT, "experiments", "colab_reap.py")
ADOPT = os.path.join(ROOT, "experiments", "colab_adopt.py")


def adopt():
    """Re-attach running VMs whose local session was lost and restart dead
    keep-alive processes (see colab_adopt.py)."""
    for acc, home in ACCOUNTS.items():
        names = [s for s, a in SESSIONS.items() if a == acc]
        try:
            r = subprocess.run([COLAB_PY, ADOPT] + names, env=dict(os.environ, HOME=home),
                               capture_output=True, text=True, timeout=180)
            for line in r.stdout.splitlines():
                if line.startswith("adopted"):
                    log(line)
        except subprocess.TimeoutExpired:
            pass
STAGE = os.environ.get("FLEET_STAGE", "/root/fleet_stage")  # cc.tgz, kapoce_src.tgz, ...
POLL = 60
PUSH_EVERY = 600

ACCOUNTS = {"a1": "/root", "a2": "/root/colab2", "a3": "/root/colab3"}
SESSIONS = {"w1": "a1", "w2": "a1", "w3": "a1", "w4": "a2", "w5": "a2", "w6": "a2",
            "w7": "a3", "w8": "a3", "w9": "a3"}
SLOTS = 1  # concurrent jobs per machine (a job itself uses both cores)
BATCH = 4  # jobs handed to a machine at once and run one after another, so that
           # machines keep working while this manager is not running

PROBE = r'''
import os, json, subprocess, glob
if (not os.path.exists("/content/SETUP_OK") and os.path.exists("/content/kapoce/build/ClusterEditing")
        and os.path.exists("/content/setup.log")
        and "SETUP_DONE" in open("/content/setup.log").read()):
    open("/content/SETUP_OK", "w").close()
out = {"setup": os.path.exists("/content/SETUP_OK"),
       "results": sorted(os.path.basename(p) for p in glob.glob("/content/results/*.json")),
       "running": [l for l in subprocess.run(["ps", "-eo", "args"], capture_output=True,
                   text=True).stdout.splitlines() if ("run_pair.py" in l or "run_seq.py" in l or "run_cert.py" in l) and "python3" in l]}
print("PROBE" + json.dumps(out))
'''


def log(msg):
    line = time.strftime("%Y-%m-%d %H:%M:%S ") + msg
    print(line, flush=True)
    with open(LOG, "a") as fh:
        fh.write(line + "\n")


def colab(session, args, timeout=300, stdin=None):
    env = dict(os.environ, HOME=ACCOUNTS[SESSIONS[session]])
    try:
        r = subprocess.run([COLAB, "--auth", "oauth2"] + args, input=stdin, env=env,
                           capture_output=True, text=True, timeout=timeout)
        return r.returncode, r.stdout + r.stderr
    except subprocess.TimeoutExpired:
        return -1, "timeout"


def load(path, default):
    try:
        with open(path) as fh:
            return json.load(fh)
    except Exception:
        return default


def save(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, indent=1)
    os.replace(tmp, path)


LAST_ERR = {}


def probe(s):
    for attempt in range(2):
        code, out = colab(s, ["exec", "-s", s], timeout=120, stdin=PROBE)
        for line in out.splitlines():
            if line.startswith("PROBE"):
                return json.loads(line[5:])
        LAST_ERR[s] = " ".join(out.split())[-160:]
        time.sleep(5)
    return None


def setup(s):
    """(Re)create the session and start the detached setup script."""
    colab(s, ["stop", "-s", s], timeout=120)
    # dead VMs keep their assignment (and count against the quota) until released
    env = dict(os.environ, HOME=ACCOUNTS[SESSIONS[s]])
    try:
        r = subprocess.run([COLAB_PY, REAP], env=env, capture_output=True, text=True,
                           timeout=180)
        if r.stdout.strip():
            log(f"{s}: {r.stdout.strip()}")
    except subprocess.TimeoutExpired:
        pass
    code, out = colab(s, ["new", "-s", s], timeout=600)
    if "READY" not in out:
        log(f"{s}: new failed: {out[-200:]}")
        return False
    for f in ["cc.tgz", "kapoce_src.tgz", "setup_full.sh", "run_pair.py", "run_seq.py", "run_cert.py"]:
        colab(s, ["upload", "-s", s, os.path.join(STAGE, f), f"/content/{f}"], timeout=900)
    launch = ('import subprocess\nsubprocess.Popen("nohup bash /content/setup_full.sh > '
              '/content/setup.log 2>&1 && touch /content/SETUP_OK &", shell=True)\nprint("ok")\n')
    colab(s, ["exec", "-s", s], timeout=120, stdin=launch)
    log(f"{s}: recreated, setup started")
    return True


CODE_FILES = ["cc.tgz", "run_pair.py", "run_seq.py", "run_cert.py"]


def refresh_code(s, st):
    """Upload the staged code again when it changed since the last upload."""
    ver = max(os.path.getmtime(os.path.join(STAGE, f)) for f in CODE_FILES)
    if st.get("code", 0) >= ver:
        return True
    for f in CODE_FILES:
        colab(s, ["upload", "-s", s, os.path.join(STAGE, f), f"/content/{f}"], timeout=900)
    code, out = colab(s, ["exec", "-s", s], timeout=180, stdin=(
        'import subprocess\nr = subprocess.run("cd /content/cc && tar xzf /content/cc.tgz", '
        'shell=True)\nprint("ok" if r.returncode == 0 else "fail")\n'))
    if "ok" in out:
        st["code"] = ver
        log(f"{s}: code refreshed")
        return True
    return False


SCRIPT = {"s": "run_seq", "c": "run_cert"}


def launch(s, jobs):
    """Start the jobs on machine s, one after another, in one detached shell."""
    # tags starting with "s": the sequential protocol (colab_run_seq.py); "c":
    # CertiFlip with a checked certificate (colab_run_cert.py)
    seq = "; ".join(f"python3 /content/{SCRIPT.get(j['tag'][0], 'run_pair')}.py "
                    f"{j['T']} {j['seed']} '{j['tag']}' "
                    f"'{j['graph']}' > /content/log_{j['id']}.txt 2>&1" for j in jobs)
    script = f"import subprocess\nsubprocess.Popen({('nohup bash -c ' + repr(seq) + ' > /dev/null 2>&1 &')!r}, shell=True)\nprint('ok')\n"
    code, out = colab(s, ["exec", "-s", s], timeout=120, stdin=script)
    return "ok" in out


def git_push():
    subprocess.run(["git", "-C", ROOT, "add", "results/colab"], capture_output=True, timeout=300)
    r = subprocess.run(["git", "-C", ROOT, "commit", "-q", "-m",
                        "Colab benchmark results (automatic sync)", "--", "results/colab"],
                       capture_output=True, timeout=300)
    if r.returncode == 0:
        for d in (2, 4, 8, 16):
            p = subprocess.run(["git", "-C", ROOT, "push", "-q", "origin", "main"],
                               capture_output=True, timeout=300)
            if p.returncode == 0:
                log("pushed results")
                return
            time.sleep(d)
            subprocess.run(["git", "-C", ROOT, "pull", "-q", "--rebase", "origin", "main"],
                           capture_output=True, timeout=300)
        log("push failed")


def merge_inbox(queue):
    """Add the job batches dropped into INBOX (by ``add``) to the queue."""
    os.makedirs(INBOX, exist_ok=True)
    ids = {j["id"] for j in queue}
    for f in sorted(os.listdir(INBOX)):
        path = os.path.join(INBOX, f)
        batch = load(path, None)
        if batch is None:
            continue
        for pre in batch.get("cancel", []):
            n0 = len(queue)
            queue[:] = [j for j in queue if not (j["status"] == "pending" and
                                                   j["id"].startswith(pre))]
            log(f"cancelled {n0 - len(queue)} pending jobs matching {pre}")
        jobs = [j for j in batch.get("jobs", []) if j["id"] not in ids]
        if batch.get("front"):
            queue[:0] = jobs
        else:
            queue.extend(jobs)
        os.remove(path)
        log(f"queued {len(jobs)} jobs from {f}")
    return queue


def step(state):
    adopt()
    queue = load(QUEUE, [])
    n0 = len(queue)
    queue = merge_inbox(queue)
    by_id = {j["id"]: j for j in queue}
    # push queue changes at once: the working copy may not survive a restart
    new_files = len(queue) != n0
    with ThreadPoolExecutor(len(SESSIONS)) as ex:
        probes = dict(zip(SESSIONS, ex.map(probe, SESSIONS)))
    if all(p is None for p in probes.values()):
        # every machine unreachable at once: either a local network problem or all
        # VMs were reclaimed (e.g. while this manager was down); ask the server
        code, out = colab(next(iter(SESSIONS)), ["sessions"], timeout=120)
        if code != 0:
            log("all probes failed and the server is unreachable; skipping this round")
            return False
        log("all probes failed but the server answers: recreating machines")
    for s in SESSIONS:
        st = state.setdefault(s, {"fails": 0, "setup_started": 0})
        pr = probes[s]
        if pr is None:
            st["fails"] += 1
            st.setdefault("fail_since", time.time())
            log(f"{s}: probe failed ({st['fails']}): {LAST_ERR.get(s, '')}")
            # a failed probe is often a stale token after this manager was down;
            # recreating destroys running jobs, so only do it after a long outage
            # "not found" after adopt(): no local session and no VM to re-attach
            gone = "not found" in LAST_ERR.get(s, "") and st["fails"] >= 2
            if gone or (st["fails"] >= 6 and time.time() - st["fail_since"] > 1800):
                for j in queue:
                    if j.get("machine") == s and j["status"] == "running":
                        j["status"], j["machine"] = "pending", None
                if setup(s):
                    st["setup_started"] = time.time()
                    st["code"] = time.time()
                st["fails"] = 0
                st.pop("fail_since", None)
            continue
        st["fails"] = 0
        st.pop("fail_since", None)
        # collect finished results
        for f in pr["results"]:
            jid = f[:-5]
            local = os.path.join(OUT, f)
            if not os.path.exists(local):
                colab(s, ["download", "-s", s, f"/content/results/{f}", local], timeout=300)
                if os.path.exists(local):
                    new_files = True
                    log(f"{s}: got {f}")
            if jid in by_id and by_id[jid]["status"] != "done":
                by_id[jid]["status"] = "done"
        if not pr["setup"]:
            # a machine that never started its setup (fresh session) gets one
            if time.time() - st.get("setup_started", 0) > 1800:
                if setup(s):
                    st["setup_started"] = time.time()
                    st["code"] = time.time()
            continue
        busy = [j for j in queue if j.get("machine") == s and j["status"] == "running"]
        # a job whose process vanished without a result is re-queued
        if busy and not pr["running"]:
            for j in busy:
                if j["id"] + ".json" not in pr["results"]:
                    log(f"{s}: job {j['id']} vanished, re-queued")
                    j["status"], j["machine"] = "pending", None
            busy = []
        if len(busy) < SLOTS and not pr["running"]:
            nxt = [j for j in queue if j["status"] == "pending"][:BATCH]
            if nxt and refresh_code(s, st) and launch(s, nxt):
                for j in nxt:
                    j["status"], j["machine"], j["started"] = "running", s, time.time()
                log(f"{s}: launched {', '.join(j['id'] for j in nxt)}")
    save(QUEUE, queue)
    save(STATE, state)
    return new_files


def main():
    os.makedirs(OUT, exist_ok=True)
    state = load(STATE, {})
    last_push = 0
    while True:
        try:
            if step(state) or time.time() - last_push > PUSH_EVERY:
                git_push()
                last_push = time.time()
        except Exception as exc:  # never die on a transient error
            log(f"error: {exc!r}")
        with open(os.path.join(OUT, "heartbeat"), "w") as fh:
            fh.write(str(time.time()))
        time.sleep(POLL)


def add(tag, T, seeds, graphs, front=False):
    jobs = [dict(id=f"{tag}_{gname.replace('/', '_')}_{seed}", tag=tag, T=T, seed=seed, graph=gname,
                 status="pending", machine=None) for seed in seeds for gname in graphs]
    os.makedirs(INBOX, exist_ok=True)
    save(os.path.join(INBOX, f"{time.time():.6f}.json"), {"front": front, "jobs": jobs})
    print(len(jobs), "jobs submitted")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "cancel":
        # cancel PREFIX...   (pending jobs whose id starts with a prefix)
        os.makedirs(INBOX, exist_ok=True)
        save(os.path.join(INBOX, f"{time.time():.6f}.json"), {"cancel": sys.argv[2:]})
    elif len(sys.argv) > 1 and sys.argv[1] in ("add", "add-front"):
        # add[-front] TAG T SEEDS(comma) GRAPH...
        add(sys.argv[2], float(sys.argv[3]), [int(x) for x in sys.argv[4].split(",")],
            sys.argv[5:], front=sys.argv[1] == "add-front")
    else:
        main()

"""Re-attach Colab VMs after a restart of this machine.

The colab CLI forgets a session when its runtime token has expired, which
happens about an hour after the keep-alive process that refreshes the token has
died (e.g. when this container was restarted), although the VM itself is still
running.  For every VM assignment of the account that no local session refers
to, this script re-creates the local session entry (fresh token and URL from the
assignment list) under the name recorded in ENDPOINTS (or a free name from
argv), and it restarts the keep-alive process of every session whose process is
gone.  Run with the Python interpreter of the colab CLI and HOME set to the
account; argv: the session names of the account.
"""
import json
import os
import sys

from colab_cli.auth import AuthProvider
from colab_cli.commands.session import spawn_keep_alive
from colab_cli.common import State
from colab_cli.state import SessionState

ENDPOINTS = os.environ.get("FLEET_ENDPOINTS", "/root/fleet_stage/endpoints.json")


def alive(pid):
    if not pid:
        return False
    try:
        os.kill(pid, 0)
        return True
    except OSError:
        return False


def main(names):
    state = State()
    try:
        emap = json.load(open(ENDPOINTS))
    except Exception:
        emap = {}
    local = state.store.list()
    known = {s.endpoint for s in local.values()}
    free = [n for n in names if n not in local]
    for a in state.client.list_assignments():
        if a.endpoint in known:
            continue
        name = emap.get(a.endpoint)
        if name not in free:
            name = next((n for n in free if n not in emap.values()), free[0] if free else None)
        if name is None:
            continue
        free.remove(name)
        s = SessionState(name=name, token=a.runtime_proxy_info.token,
                         url=a.runtime_proxy_info.url, endpoint=a.endpoint)
        state.store.add(s)
        emap[a.endpoint] = name
        print("adopted", name, a.endpoint)
    for name, s in state.store.list().items():
        emap[s.endpoint] = name
        if not alive(s.keep_alive_pid):
            s.keep_alive_pid = spawn_keep_alive(s.endpoint, name, auth_provider=AuthProvider.OAUTH2)
            state.store.add(s)
            print("keep-alive", name)
    tmp = ENDPOINTS + ".tmp"
    json.dump(emap, open(tmp, "w"), indent=1)
    os.replace(tmp, ENDPOINTS)


if __name__ == "__main__":
    main(sys.argv[1:])

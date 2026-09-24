"""Release Colab VM assignments that no local session refers to (VMs whose
kernel died but whose assignment still counts against the session quota).
Run with the Python interpreter of the colab CLI, with HOME set to the account."""
from colab_cli.common import State

state = State()
known = {s.endpoint for s in (state.store.get(n) for n in state.store.list_names())} \
    if hasattr(state.store, "list_names") else None
if known is None:
    import json, os
    path = os.path.expanduser("~/.config/colab-cli/sessions.json")
    known = {v["endpoint"] for v in json.load(open(path)).values()} if os.path.exists(path) else set()
for a in state.client.list_assignments():
    if a.endpoint not in known:
        state.client.unassign(a.endpoint)
        print("released", a.endpoint)

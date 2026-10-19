import sys, itertools, ast
from w6a_exh import *
pats4 = [ast.literal_eval(l) for l in open(sys.argv[1]).read().strip().split("\n")]
m = int(sys.argv[2]); kx = int(sys.argv[3]); part = int(sys.argv[4]); nparts = int(sys.argv[5])
cands = [p + (q,) for p in pats4 for q in itertools.product('ec', 'ec')]
for i, pat in enumerate(cands):
    if i % nparts != part: continue
    res = run(m, kx, pat)
    print(m, kx, pat, len(res), flush=True)
    for r in res[:2]: print("   ", r, flush=True)

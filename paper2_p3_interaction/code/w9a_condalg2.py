import sys, random, itertools, time
sys.path.insert(0,'/home/user/correlation_clustering/paper2_p3_interaction/code')
from w9_lib import *
from w9a_condalg import *
rng=random.Random(int(sys.argv[1])); N=int(sys.argv[2]); l=int(sys.argv[3])
cnt=0; mism=0; hist={}; nodes=0; needfar=0
while cnt<N:
    m=rng.randint(3,4)
    n0,E0,soft,H=gen_path(rng,m,extra=rng.randint(0,2),topo=rng.choice(['path','tree']))
    if n0>13: continue
    hardpairs=[key(u,v) for u,v in itertools.combinations(range(n0),2) if key(u,v) not in soft]
    # random sparse hard edges (chords) kept statically feasible
    cur=set()
    for _ in range(rng.randint(1,5)):
        p=rng.choice(hardpairs); new=set(cur)
        new^={p}
        inst=Inst(n0,E0|new,H)
        if not any(len(c[0])==0 for c in inst.compile()): cur=new
    inst=Inst(n0,E0|cur,H)
    if not inst.valid_packing(): continue
    A=Alg(inst)
    O,_=exact_opt(inst.n,sorted(inst.E))
    ex=O-len(H)
    if ex==0: continue
    cnt+=1
    ans=A.solve(frozenset(),l)
    nodes+=A.nodes
    truth=(ex<=l)
    # does the instance need far edits? opt with U=empty (breaks+flips only) > l but truth True
    v0=A.opt(frozenset(),l)
    if truth and v0>l: needfar+=1
    hist[(ex,ans)]=hist.get((ex,ans),0)+1
    if ans!=truth:
        mism+=1; print('MISMATCH',ex,ans,inst.n,sorted(inst.E),H,flush=True)
print('instances',cnt,'mismatches',mism,'needing far edits',needfar,'hist',sorted(hist.items()),'avg nodes',nodes/cnt)

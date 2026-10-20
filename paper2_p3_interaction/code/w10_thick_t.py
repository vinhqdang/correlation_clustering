import sys, time, random, itertools
from w10_red import *
k,n,seed,dens,p=int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]),float(sys.argv[4]),int(sys.argv[5])
rng=random.Random(seed)
S=[[set() for _ in range(k)] for _ in range(k)]
for i in range(k):
    for j in range(k):
        for a in range(1,n+1):
            for b in range(1,n+1):
                if rng.random()<dens: S[i][j].add((a,b))
        if not S[i][j]: S[i][j].add((1,1))
t0=time.time()
W,meta=build_gt(k,n,S,p,thin_mode=False)
I=Inst(W)
print('thick n',I.n,'|H|',len(I.members),I.valid,'c',I.c,'max matchings used',max(W.used.values()),'of',p-1,round(time.time()-t0,1),flush=True)
sols=gt_solvable(k,n,S)
cells=[(i,j) for i in range(k) for j in range(k)]
choices=list(itertools.product(*[sorted(S[i][j]) for (i,j) in cells]))
rng.shuffle(choices)
cnt={True:0,False:0}
for choice in choices:
    t={c:e for c,e in zip(cells,choice)}
    exp=t in sols
    if cnt[exp]>=int(sys.argv[6]): continue
    cnt[exp]+=1
    t1=time.time()
    lab=I.solve(drop_hard=cut_choice(meta,t))
    print('choice',choice,'GT-valid',exp,'ZE feasible',lab is not None,round(time.time()-t1,1),flush=True)

from w8a_lib import *
import collections
for m in (3,4,5):
    n=2*m+1
    E=[(i,i+1) for i in range(2*m)]+[(1,2*m-1)]
    g=G8(n,E)
    opt,sols=g.opt()
    parts=[part_from_mask(s,n) for s in sols]
    edl=[edits_of_partition(g,p) for p in parts]
    packs=enum_packings(g,minsize=m)
    print('FAR m=%d  optima=%d  packings of size>=m: %d'%(m,len(parts),len(packs)))
    stat=collections.Counter()
    best=None
    for H in packs:
        H=list(H); pm=pair_maximal(g,H)
        I=Igraph(H); tw=tw_exact(I)
        # best optimum: min over optima of max over hard edits of distance (inf for disconnected)
        rbest=None
        for S in edl:
            he=hard_edits(H,S)
            r=0
            for (a,b) in he:
                d=member_dist(H,a,b,I)
                r=max(r, 0 if d==-1 else d)
            if rbest is None or r<rbest: rbest=r
        stat[(len(H),pm,tw,rbest)]+=1
        if len(H)==m and pm and rbest is not None and (best is None or rbest<best[0]): best=(rbest,H)
    for k,v in sorted(stat.items()): print('  (|H|,pairmax,tw(I),min-over-optima max hard-edit distance) =',k,'count',v)
    print('  best repack:',best)

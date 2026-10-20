"""exact test M1: k=1 reduction (chains A,B + tunnels), Grid Tiling k=1,n=1 (extra end pad to reach m>=8), p given; l=0,1."""
import sys, time
from w10_red import *
from w10_exact import exact_decide
p=int(sys.argv[1]); mode=sys.argv[2]
n=1
U,m,sg,sep=layout2(n,1)
print('m',m,'U',U, 'slot group', sg)
tilesA={sg[(1,1)]} if mode=='yes' else set()
tilesB={sg[(1,1)]} if mode=='yes' else set()
tunnels=[((('A',0,0),sg[(1,1)]),(('B',0,0),sg[(1,1)]))] if mode=='yes' else []
W,meta=build_general(n,p,{('A',0,0):tilesA,('B',0,0):tilesB},[],tunnels,extra=1)
I=Inst(W)
print('thick n',I.n,'|H|',len(I.members),I.valid,'c',I.c,'matchings',max(W.used.values()),flush=True)
for l in (0,1,2):
    r,sol=exact_decide(I.n,I.E,I.members,l,timeout=3000)
    print('l',l,'OPT<=|H|+l ?',r,flush=True)

"""exact test M2: row gadget = two A-chains (n=2) with a comparator, each anchor tunnelled to a private stub; budget l=2.
YES iff exists tiles e1 in T1, e2 in T2 with a(e1) <= a(e2).   usage: w10_ex2.py p "T1" "T2" ltest   e.g.  "[(1,1)]" "[(2,1)]" """
import sys, time, ast
from w10_red import *
from w10_exact import exact_decide
p=int(sys.argv[1]); T1=ast.literal_eval(sys.argv[2]); T2=ast.literal_eval(sys.argv[3]); ls=[int(x) for x in sys.argv[4].split(',')]
n=2
U,m,sg,sep=layout2(n)
tiles={'A1':{sg[e] for e in T1}, 'A2':{sg[e] for e in T2}}
tun=[(('A1',sg[e]),None) for e in T1]+[(('A2',sg[e]),None) for e in T2]
W,meta=build_general(n,p,tiles,[('A1','A2',1)],tun)
I=Inst(W)
exp=any(e1[0]<=e2[0] for e1 in T1 for e2 in T2)
print('T1',T1,'T2',T2,'expected YES(l=2):',exp,'| thick n',I.n,'|H|',len(I.members),I.valid,'c',I.c,flush=True)
for l in ls:
    t0=time.time()
    r,sol=exact_decide(I.n,I.E,I.members,l,timeout=20000)
    print('RESULT p',p,'T1',T1,'T2',T2,'l',l,'OPT<=|H|+l ?',r,round(time.time()-t0,1),flush=True)

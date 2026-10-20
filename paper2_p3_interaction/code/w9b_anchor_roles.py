import pickle, re, ast
from w9b_k3search import chain_inst
from w9b_pend import LIST2ROLE
good=pickle.load(open('/tmp/w9b_k3good.pkl','rb'))
fam_by_idx={}
for line in open('/tmp/k3fam.log'):
    m=re.match(r'(\d+) sig .*family \(member,list\): (.*)$',line.strip())
    if m: fam_by_idx[int(m.group(1))]=ast.literal_eval(m.group(2))
for idx,(R,sig,cfg,outm,inm,eps,ns,k2) in enumerate(good):
    if idx not in fam_by_idx: continue
    I=chain_inst(4,cfg,outm,inm,eps)
    out=[]
    for m,L in fam_by_idx[idx]:
        if len(L)!=2: continue
        role,val=LIST2ROLE[tuple(sorted(L))]
        a,b,c=I.H[m]
        pr={0:(a,b),1:(b,c),2:(a,c)}[role]
        info=[]
        for x in pr:
            roles=[]
            for mm in I.mem[x]:
                aa,bb,cc=I.H[mm]
                roles.append('centre' if bb==x else 'end')
            info.append(roles)
        out.append((m,pr,info))
        break
    print(idx,'sig',sig,'guard member',out[0][0],'guarded pair',out[0][1],'anchor roles',out[0][2])

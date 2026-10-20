from w8a_lib import *
for m in (3,4,5,6):
    n=2*m+1
    E=[(i,i+1) for i in range(2*m)]+[(1,2*m-1)]
    g=G8(n,E)
    nu,Hn=nu_mod(g)
    opt,_=g.opt(enum=False)
    print('FAR m=%d n=%d |E|=%d OPT=%d nu_mod=%d  (nu<=floor(|E|/2)=%d)'%(m,n,len(E),opt,nu,len(E)//2))

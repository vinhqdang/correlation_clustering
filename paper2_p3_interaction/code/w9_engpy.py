import ctypes, os, numpy as np
_lib=ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)),'w9_eng.so'))
_lib.eval.argtypes=[ctypes.c_int,ctypes.c_char_p,ctypes.c_int,ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_int)]
_lib.eval.restype=ctypes.c_int
_lib.get_e.argtypes=[ctypes.c_int]; _lib.get_e.restype=ctypes.c_int
def evaluate(n,E,H):
    adj=bytearray(n*n)
    for u,v in E:
        adj[u*n+v]=1; adj[v*n+u]=1
    Hin=(ctypes.c_int*(3*len(H)))(*[x for t in H for x in t])
    out=(ctypes.c_int*4)()
    rc=_lib.eval(n,bytes(adj),len(H),Hin,out)
    if rc!=0: return rc,None
    return 0,(out[0],out[1])
def e_table(m):
    return [None]+[_lib.get_e(R) for R in range(1,1<<m)]
_lib.get_er.argtypes=[ctypes.c_int,ctypes.c_int]; _lib.get_er.restype=ctypes.c_int
def e_relaxed(R,relax): return _lib.get_er(R,relax)
_lib.build_only.argtypes=[ctypes.c_int,ctypes.c_char_p,ctypes.c_int,ctypes.POINTER(ctypes.c_int)]
_lib.build_only.restype=ctypes.c_int
_lib.single_relax.argtypes=[ctypes.POINTER(ctypes.c_int)]; _lib.single_relax.restype=ctypes.c_int
_lib.cost_full_relaxed.argtypes=[ctypes.c_int]; _lib.cost_full_relaxed.restype=ctypes.c_int
def single_relax(n,E,H):
    adj=bytearray(n*n)
    for u,v in E:
        adj[u*n+v]=1; adj[v*n+u]=1
    Hin=(ctypes.c_int*(3*len(H)))(*[x for t in H for x in t])
    rc=_lib.build_only(n,bytes(adj),len(H),Hin)
    if rc!=0: return rc,None,None
    out=(ctypes.c_int*len(H))()
    base=_lib.single_relax(out)
    return 0,base,list(out)
def cost_relaxed(Z):
    return _lib.cost_full_relaxed(Z)

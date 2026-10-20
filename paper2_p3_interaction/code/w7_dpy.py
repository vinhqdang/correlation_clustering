import ctypes, os, numpy as np
_lib=ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)),'w7_dp.so'))
_lib.solve.argtypes=[ctypes.c_int,ctypes.POINTER(ctypes.c_uint),ctypes.POINTER(ctypes.c_int),ctypes.c_int,ctypes.POINTER(ctypes.c_int)]
_lib.solve.restype=ctypes.c_int
def opt(n,edges,maxsol=0):
    adj=(ctypes.c_uint*n)()
    for u,v in edges: adj[u]|=1<<v; adj[v]|=1<<u
    lab=(ctypes.c_int*(max(1,maxsol)*n))(); ns=ctypes.c_int(0)
    r=_lib.solve(n,adj,lab,maxsol,ctypes.byref(ns))
    sols=[tuple(lab[i*n:(i+1)*n]) for i in range(ns.value)]
    return r,sols

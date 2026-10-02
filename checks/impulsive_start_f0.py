# Part of the validation of the lee-wave benchmark. Own calculation, not a published result.
# Needs numpy and scipy:  pip install numpy scipy
# Linear, inviscid, f=0 NONHYDROSTATIC response to an IMPULSIVE START of flow U over the Gaussian ridge
# (tall domain with rigid lid at H so that no wave returns within 6 h).
# Per wavenumber k, in the ridge frame (Baker & Mashayek's steady limit is the t->infinity answer):
#   q = psi_zz - k^2 psi ;  q_t = -i k U q + i k b ;  b_t = -i k U b - i k N^2 psi
#   psi(0) = U*hhat(k) (bottom BC w = U h_x, switched on at t=0), psi(H)=0.   u=-psi_z, w = i k psi.
import numpy as np
from scipy.linalg import solve_banded
N,U,h0,sig=2e-3,0.2,20.,500.
H,nz=3000.,750; dz=H/nz; z=np.arange(nz+1)*dz
L=200e3; dk=2*np.pi/L; ks=np.arange(1,int(0.008/dk)+1)*dk
hhat=h0*sig*np.sqrt(2*np.pi)*np.exp(-ks**2*sig**2/2)
ref=N*U*h0**2
dt=100.; T=6*3600.; nsteps=int(T/dt)
def make_solver(k):
    n=nz-1; ab=np.zeros((3,n))
    ab[0,1:]=1/dz**2; ab[1,:]=-2/dz**2-k*k; ab[2,:-1]=1/dz**2
    def solve(q,psi0):
        rhs=q[1:-1].copy(); rhs[0]-=psi0/dz**2     # psi(0)=psi0, psi(H)=0
        psi=np.zeros(nz+1,complex); psi[0]=psi0
        psi[1:-1]=solve_banded((1,1),ab,rhs); return psi
    return solve
solvers=[make_solver(k) for k in ks]
hts=[25,50,75,100,150,200,300,400,500]
jm=[int(round((20.+hh)/dz)) for hh in hts]
Q=np.zeros((len(ks),nz+1),complex); B=np.zeros_like(Q)
def rhs(Q,B):
    dQ=np.empty_like(Q); dB=np.empty_like(B); PS=np.empty_like(Q)
    for i,k in enumerate(ks):
        psi=solvers[i](Q[i],U*hhat[i]); PS[i]=psi
        dQ[i]=-1j*k*U*Q[i]+1j*k*B[i]
        dB[i]=-1j*k*U*B[i]-1j*k*N**2*psi
    return dQ,dB,PS
def flux(PS):
    out=[]
    for j in jm:
        u=-np.gradient(PS[:,:],dz,axis=1)[:,j]; w=1j*ks*PS[:,j]
        out.append(-np.sum(2*np.real(u*np.conj(w)))*dk/(2*np.pi)/ref)
    return out
rec=[]
for n in range(nsteps+1):
    t=n*dt
    if n%int(900/dt)==0:
        _,_,PS=rhs(Q,B); rec.append((t/3600,*flux(PS)))
    k1q,k1b,_=rhs(Q,B)
    k2q,k2b,_=rhs(Q+0.5*dt*k1q,B+0.5*dt*k1b)
    k3q,k3b,_=rhs(Q+0.5*dt*k2q,B+0.5*dt*k2b)
    k4q,k4b,_=rhs(Q+dt*k3q,B+dt*k3b)
    Q=Q+dt/6*(k1q+2*k2q+2*k3q+k4q); B=B+dt/6*(k1b+2*k2b+2*k3b+k4b)

sel=[r for r in rec if 4.0-1e-9<=r[0]<=6.0+1e-9]
avg=np.mean([r[1:] for r in sel],axis=0)
print("Linear (inviscid, nonhydrostatic, f=0) impulsive-start response, mean over t=4-6 h; steady unbounded value = 0.980")
print("height above crest [m]:", hts)
print("drag/ref               :", np.round(avg,3).tolist())

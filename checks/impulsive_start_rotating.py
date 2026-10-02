# Part of the validation of the lee-wave benchmark. Own calculation, not a published result.
# Needs numpy and scipy:  pip install numpy scipy
# Linear, inviscid, NONHYDROSTATIC impulsive-start response WITH rotation (f-plane), tall rigid-lid domain.
# Per wavenumber k (ridge frame, Dt = d/dt + i k U):
#   q = psi_zz - k^2 psi ;  Dt q = i k b - f v_z ;  Dt v = f psi_z ;  Dt b = -i k N^2 psi
#   psi(0) = U hhat(k) ; psi(H) = 0 ; u = -psi_z, w = i k psi.   Conserved flux F = <u w> - f <v b>/N^2.
import numpy as np, sys
from scipy.linalg import solve_banded
N,U,h0,sig=2e-3,0.2,20.,500.
H,nz=3000.,750; dz=H/nz
L=200e3; dk=2*np.pi/L; ks=np.arange(1,int(0.008/dk)+1)*dk
hhat=h0*sig*np.sqrt(2*np.pi)*np.exp(-ks**2*sig**2/2); ref=N*U*h0**2
hts=[25,50,75,100,150,200,300]; jm=[int(round((20.+h)/dz)) for h in hts]
dt=100.; T=6*3600.; nsteps=int(T/dt)
def make_solver(k):
    n=nz-1; ab=np.zeros((3,n)); ab[0,1:]=1/dz**2; ab[1,:]=-2/dz**2-k*k; ab[2,:-1]=1/dz**2
    def solve(q,psi0):
        rhs=q[1:-1].copy(); rhs[0]-=psi0/dz**2
        psi=np.zeros(nz+1,complex); psi[0]=psi0; psi[1:-1]=solve_banded((1,1),ab,rhs); return psi
    return solve
solvers=[make_solver(k) for k in ks]
def run(f):
    Q=np.zeros((len(ks),nz+1),complex); V=np.zeros_like(Q); B=np.zeros_like(Q)
    def rhs(Q,V,B):
        dQ=np.empty_like(Q); dV=np.empty_like(V); dB=np.empty_like(B); PS=np.empty_like(Q)
        for i,k in enumerate(ks):
            psi=solvers[i](Q[i],U*hhat[i]); PS[i]=psi
            psiz=np.gradient(psi,dz); vz=np.gradient(V[i],dz)
            dQ[i]=-1j*k*U*Q[i]+1j*k*B[i]-f*vz
            dV[i]=-1j*k*U*V[i]+f*psiz
            dB[i]=-1j*k*U*B[i]-1j*k*N**2*psi
        return dQ,dV,dB,PS
    def drag(PS,V,B):
        out=[]
        P=lambda a,c: np.sum(2*np.real(a*np.conj(c)))*dk/(2*np.pi)
        for j in jm:
            u=-np.gradient(PS,dz,axis=1)[:,j]; w=1j*ks*PS[:,j]
            uw=P(u,w); vb=P(V[:,j],B[:,j]); out.append(-(uw-f*vb/N**2)/ref)
        return out
    rec=[]
    for n in range(nsteps+1):
        t=n*dt
        if n%int(900/dt)==0 and t>=4*3600-1:
            _,_,_,PS=rhs(Q,V,B); rec.append(drag(PS,V,B))
        a=rhs(Q,V,B); b=rhs(Q+.5*dt*a[0],V+.5*dt*a[1],B+.5*dt*a[2])
        c=rhs(Q+.5*dt*b[0],V+.5*dt*b[1],B+.5*dt*b[2]); d=rhs(Q+dt*c[0],V+dt*c[1],B+dt*c[2])
        Q=Q+dt/6*(a[0]+2*b[0]+2*c[0]+d[0]); V=V+dt/6*(a[1]+2*b[1]+2*c[1]+d[1]); B=B+dt/6*(a[2]+2*b[2]+2*c[2]+d[2])
    return np.mean(rec,axis=0)
d0=run(0.0); d1=run(1.4e-4)
print("height above crest [m]      :",hts)
print("linear transient f=0      /ref:",np.round(d0,3).tolist())
print("linear transient f=1.4e-4 /ref:",np.round(d1,3).tolist())
print("ratio rot/nonrot (linear)     :",np.round(d1/d0,3).tolist())

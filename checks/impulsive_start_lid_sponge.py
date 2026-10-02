# Part of the validation of the lee-wave benchmark. Own calculation, not a published result.
# Needs numpy and scipy:  pip install numpy scipy
# Linear impulsive-start response in the SIMULATION'S OWN DOMAIN: rigid lid at H = 1000 m and the Rayleigh sponge
# (top 200 m, 600 s, quadratic) damping u, w and b'.  f = 0, nonhydrostatic.  (x-sponges ignored.)
#   q = psi_zz - k^2 psi ;  Dt q = i k b - r q - r_z psi_z ;  Dt b = -i k N^2 psi - r b ;  Dt = d/dt + i k U
import numpy as np
from scipy.linalg import solve_banded
N,U,h0,sig=2e-3,0.2,20.,500.
H,nz=1000.,500; dz=H/nz; z=np.arange(nz+1)*dz
def run(top=200.,tau=600.):
    r=np.where(z>H-top,((z-(H-top))/top)**2,0.0)/tau; rz=np.gradient(r,dz)
    L=200e3; dk=2*np.pi/L; ks=np.arange(1,int(0.008/dk)+1)*dk
    hhat=h0*sig*np.sqrt(2*np.pi)*np.exp(-ks**2*sig**2/2); ref=N*U*h0**2
    hts=[25,50,75,100,150,200,300,400,500]; jm=[int(round((20.+h)/dz)) for h in hts]
    dt=100.; nsteps=int(6*3600/dt)
    def mk(k):
        n=nz-1; ab=np.zeros((3,n)); ab[0,1:]=1/dz**2; ab[1,:]=-2/dz**2-k*k; ab[2,:-1]=1/dz**2
        def solve(q,p0):
            rhs=q[1:-1].copy(); rhs[0]-=p0/dz**2
            psi=np.zeros(nz+1,complex); psi[0]=p0; psi[1:-1]=solve_banded((1,1),ab,rhs); return psi
        return solve
    S=[mk(k) for k in ks]
    Q=np.zeros((len(ks),nz+1),complex); B=np.zeros_like(Q)
    def rhs(Q,B):
        dQ=np.empty_like(Q); dB=np.empty_like(B); PS=np.empty_like(Q)
        for i,k in enumerate(ks):
            psi=S[i](Q[i],U*hhat[i]); PS[i]=psi; psiz=np.gradient(psi,dz)
            dQ[i]=-1j*k*U*Q[i]+1j*k*B[i]-r*Q[i]-rz*psiz
            dB[i]=-1j*k*U*B[i]-1j*k*N**2*psi-r*B[i]
        return dQ,dB,PS
    def drag(PS):
        out=[]
        for j in jm:
            u=-np.gradient(PS,dz,axis=1)[:,j]; w=1j*ks*PS[:,j]
            out.append(-np.sum(2*np.real(u*np.conj(w)))*dk/(2*np.pi)/ref)
        return out
    rec=[]; series=[]
    for n in range(nsteps+1):
        t=n*dt
        if n%int(900/dt)==0:
            _,_,PS=rhs(Q,B); d=drag(PS); series.append((t/3600,d[2]))
            if t>=4*3600-1: rec.append(d)
        a=rhs(Q,B); b=rhs(Q+.5*dt*a[0],B+.5*dt*a[1]); c=rhs(Q+.5*dt*b[0],B+.5*dt*b[1]); d_=rhs(Q+dt*c[0],B+dt*c[1])
        Q=Q+dt/6*(a[0]+2*b[0]+2*c[0]+d_[0]); B=B+dt/6*(a[1]+2*b[1]+2*c[1]+d_[1])
    return hts,np.mean(rec,axis=0),series
hts,avg,series=run()
print("Linear impulsive start, H=1000 m rigid lid + simulation sponge, f=0, mean over t=4-6 h")
print("height above crest [m]:",hts)
print("drag/ref              :",np.round(avg,3).tolist())
print("time series at 75 m (t[h], drag/ref):",[(round(a,2),round(b,3)) for a,b in series[::4]])

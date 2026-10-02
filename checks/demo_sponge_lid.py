# Part of the validation of the lee-wave benchmark. Own calculation, not a published result.
# Needs numpy and scipy:  pip install numpy scipy
# Does the DEMO's sponge (damping u, v, w but NOT buoyancy) plus a rigid lid explain a flux well above linear theory?
# Linear, f=0, nonhydrostatic, H=1000 m rigid lid, sponge in the top 200 m (tau=600 s, quadratic). Flux at z=-600 m (400 m above the floor).
import numpy as np
from scipy.linalg import solve_banded
N,U,h0,sig,H=2e-3,0.2,20.,500.,1000.
ref=N*U*h0**2
def hhat(k): return h0*sig*np.sqrt(2*np.pi)*np.exp(-k**2*sig**2/2)
# ---------- steady state, damping on u,w (+ optionally b)
def steady(damp_b, top=200., tau=600., nz=2000, Lbox=400e3, kmax=0.012):
    dz=H/nz; z=np.arange(nz+1)*dz; zh=0.5*(z[1:]+z[:-1])
    r=lambda zz: np.where(zz>H-top, ((zz-(H-top))/top)**2, 0.0)/tau
    dk=2*np.pi/Lbox; ks=np.arange(1,int(kmax/dk)+1)*dk; tot=0.0
    j=int(round(400./dz))
    for k in ks:
        s=1j*k*U+r(z); sh=1j*k*U+r(zh)
        bterm=(N**2/s) if damp_b else (N**2/(1j*k*U))*np.ones_like(s)
        a=sh[:-1]/(k*k*dz*dz); c=sh[1:]/(k*k*dz*dz); b=-(a+c)-(s[1:-1]+bterm[1:-1])
        w0=1j*k*U*hhat(k); rhs=np.zeros(nz-1,complex); rhs[0]=-a[0]*w0
        ab=np.zeros((3,nz-1),complex); ab[0,1:]=c[:-1]; ab[1,:]=b; ab[2,:-1]=a[1:]
        w=np.zeros(nz+1,complex); w[0]=w0; w[1:-1]=solve_banded((1,1),ab,rhs)
        u=1j*np.gradient(w,dz)/k
        tot+=2*np.real(u[j]*np.conj(w[j]))*dk/(2*np.pi)
    return -tot/ref
print("STEADY drag/ref at z=-600 m:  sponge on u,w,b (as in run_case.jl): %.3f   |   sponge on u,w only (as in the demo): %.3f"%(steady(True),steady(False)))
# ---------- transient, 18 h, sponge on u,w only (demo)
nz=500; dz=H/nz; z=np.arange(nz+1)*dz
top,tau=200.,600.
r=np.where(z>H-top,((z-(H-top))/top)**2,0.0)/tau; rz=np.gradient(r,dz)
L=200e3; dk=2*np.pi/L; ks=np.arange(1,int(0.008/dk)+1)*dk
hh=np.array([hhat(k) for k in ks]); dt=100.; nsteps=int(18*3600/dt); j=int(round(400./dz))
def mk(k):
    n=nz-1; ab=np.zeros((3,n)); ab[0,1:]=1/dz**2; ab[1,:]=-2/dz**2-k*k; ab[2,:-1]=1/dz**2
    def solve(q,p0):
        rhs=q[1:-1].copy(); rhs[0]-=p0/dz**2
        psi=np.zeros(nz+1,complex); psi[0]=p0; psi[1:-1]=solve_banded((1,1),ab,rhs); return psi
    return solve
S=[mk(k) for k in ks]
Q=np.zeros((len(ks),nz+1),complex); B=np.zeros_like(Q)
def rhs(Q,B,damp_b):
    dQ=np.empty_like(Q); dB=np.empty_like(B); PS=np.empty_like(Q)
    for i,k in enumerate(ks):
        psi=S[i](Q[i],U*hh[i]); PS[i]=psi; pz=np.gradient(psi,dz)
        dQ[i]=-1j*k*U*Q[i]+1j*k*B[i]-r*Q[i]-rz*pz
        dB[i]=-1j*k*U*B[i]-1j*k*N**2*psi-(r*B[i] if damp_b else 0)
    return dQ,dB,PS
def flux(PS):
    u=-np.gradient(PS,dz,axis=1)[:,j]; w=1j*ks*PS[:,j]
    return -np.sum(2*np.real(u*np.conj(w)))*dk/(2*np.pi)/ref
def transient(damp_b):
    Q=np.zeros((len(ks),nz+1),complex); B=np.zeros_like(Q); out=[]
    for n in range(nsteps+1):
        if n%int(1800/dt)==0: out.append((n*dt/3600, flux(rhs(Q,B,damp_b)[2])))
        a=rhs(Q,B,damp_b); b=rhs(Q+.5*dt*a[0],B+.5*dt*a[1],damp_b); c=rhs(Q+.5*dt*b[0],B+.5*dt*b[1],damp_b); d=rhs(Q+dt*c[0],B+dt*c[1],damp_b)
        Q=Q+dt/6*(a[0]+2*b[0]+2*c[0]+d[0]); B=B+dt/6*(a[1]+2*b[1]+2*c[1]+d[1])
    return out
o=transient(False)
print("TRANSIENT, sponge on u,w only: drag/ref at z=-600 m (t[h]: value)")
print(", ".join(f"{t:.0f}h: {v:.2f}" for t,v in o[::4]))
tt=np.array([t for t,_ in o]); vv=np.array([v for _,v in o])
print("mean over last 2 h: %.2f   (min %.2f, max %.2f)"%(vv[tt>=16].mean(),vv[tt>=16].min(),vv[tt>=16].max()))

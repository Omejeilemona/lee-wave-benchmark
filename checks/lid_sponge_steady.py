# Part of the validation of the lee-wave benchmark. Own calculation, not a published result.
# Needs numpy and scipy:  pip install numpy scipy
# Steady linear lee-wave solution WITH a rigid lid at z=H and a Rayleigh sponge in the top layer,
# to test whether lid/sponge reflection could explain the simulated drag (+48 % vs unbounded theory).
# Model: f=0, nonhydrostatic, uniform N and U, Rayleigh damping r(z) on u, w and b'.
# For wavenumber k, s(z)=i k U + r(z):  (s w_z)_z / k^2 = (s + N^2/s) w ,  w(0)=i k U hhat(k), w(H)=0
import numpy as np
from scipy.linalg import solve_banded
N,U,h0,sig,H=2e-3,0.2,20.,500.,1000.
ref=N*U*h0**2
def drag(top_thick=200., tau=600., nz=2000, Lbox=400e3, kmax=0.012, zmeas=(70.,120.), power=2, lid=True):
    dz=H/nz; z=np.arange(nz+1)*dz
    zh=0.5*(z[1:]+z[:-1])
    r =lambda zz: np.where(zz>H-top_thick, ((zz-(H-top_thick))/top_thick)**power, 0.0)/tau
    dk=2*np.pi/Lbox; ks=np.arange(1,int(kmax/dk)+1)*dk
    hhat=lambda k: h0*sig*np.sqrt(2*np.pi)*np.exp(-k**2*sig**2/2)
    tot=np.zeros(len(zmeas))
    jm=[int(round(zz/dz)) for zz in zmeas]
    for k in ks:
        s=1j*k*U+r(z); sh=1j*k*U+r(zh)
        a=sh[:-1]/(k*k*dz*dz)        # coupling to j-1 (interior j=1..nz-1)
        c=sh[1:]/(k*k*dz*dz)         # coupling to j+1
        b=-(a+c)-(s[1:-1]+N**2/s[1:-1])
        w0=1j*k*U*hhat(k)
        rhs=np.zeros(nz-1,complex); rhs[0]=-a[0]*w0
        ab=np.zeros((3,nz-1),complex); ab[0,1:]=c[:-1]; ab[1,:]=b; ab[2,:-1]=a[1:]
        w=np.zeros(nz+1,complex); w[0]=w0
        w[1:-1]=solve_banded((1,1),ab,rhs)
        if not lid:
            raise NotImplementedError
        wz=np.gradient(w,dz); u=1j*wz/k
        for i,j in enumerate(jm):
            tot[i]+=2*np.real(u[j]*np.conj(w[j]))*dk/(2*np.pi)
    return -tot/ref   # drag/ref at the measurement heights
for label,kw in [("sim sponge: 200 m, tau 600 s, quadratic",dict()),
                 ("thicker sponge 600 m, tau 600 s",dict(top_thick=600.)),
                 ("thick + strong 600 m, tau 200 s",dict(top_thick=600.,tau=200.)),
                 ("no sponge (lid only), weak damping check: 200 m, tau 1e9",dict(tau=1e9))]:
    d=drag(**kw); print(f"{label:62s} drag/ref at 70 and 120 m above floor: {d[0]:.3f}, {d[1]:.3f}")

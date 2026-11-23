#!/usr/bin/env python3
"""Fixed-q Higgs-pair shift: exact boundary and explicit contour check."""
import argparse
import hashlib
import json
import platform
from pathlib import Path
import numpy as np
import sympy as sp

ETA=np.diag([1.,-1.,-1.,-1.])
BASIS=np.eye(4)

def dot(x,y): return x@ETA@y

def action_vertex(legs, a, m, M):
    """Exp(iS) coefficient by polarization of degree-n Fourier L_n.
    A vector leg carries its full contravariant polarization.
    Scalar p is incoming even for an internal leg.
    """
    n=len(legs)
    if n not in (3,4): raise ValueError(n)
    if np.max(np.abs(sum((p for (_,p,_) in legs), start=np.zeros(4)))) > 1e-10:
        raise ValueError('vertex momentum mismatch')
    ans=0j
    v=m/a
    for mask in range(1,1<<n):
        W=np.zeros(4,dtype=complex); h=0.; chi=0.; dh=np.zeros(4,dtype=complex); dc=dh.copy()
        for i,(kind,p,e) in enumerate(legs):
            if not mask & (1<<i): continue
            if kind=='W': W+=e
            elif kind=='h': h+=1; dh-=1j*p
            elif kind=='c': chi+=1; dc-=1j*p
            else: raise ValueError(kind)
        if n==3:
            L=a*dot(W,h*dc-chi*dh)+a*m*dot(W,W)*h-M*M/(2*v)*h*(h*h+chi*chi)
        else:
            L=(a*a/2)*dot(W,W)*(h*h+chi*chi)-M*M/(8*v*v)*(h*h+chi*chi)**2
        ans+=(-1)**(n-mask.bit_count())*L
    return 1j*ans

def vector_prop(Q,m,xi):
    q2=dot(Q,Q)
    return -1j*(ETA-(1-xi)*np.outer(Q,Q)/(q2-xi*m*m))/(q2-m*m)

def direct_amplitudes(q,K,p,r,e,a,m,M,xi):
    """Returns q-contracted vector amplitude and unphased Goldstone amplitude."""
    root=[('W',q,q),('c',q,None)]
    Kleg=('W',K,e)
    out=np.zeros((2,6), dtype=complex)
    for i,off in enumerate(root):
        for j,(x,y) in enumerate(((p,r),(r,p))):
            Q=q+x
            xh=('h',x,None); yh=('h',y,None)
            left=np.array([action_vertex([off,xh,('W',-Q,u)],a,m,M) for u in BASIS])
            right=np.array([action_vertex([Kleg,yh,('W',Q,u)],a,m,M) for u in BASIS])
            out[i,2*j]=left@vector_prop(Q,m,xi)@right
            out[i,2*j+1]=(action_vertex([off,xh,('c',-Q,None)],a,m,M)*
                           (1j/(dot(Q,Q)-xi*m*m))*
                           action_vertex([Kleg,yh,('c',Q,None)],a,m,M))
        H=p+r
        out[i,4]=(action_vertex([off,Kleg,('h',H,None)],a,m,M)*
                    (1j/(dot(H,H)-M*M))*
                    action_vertex([('h',p,None),('h',r,None),('h',-H,None)],a,m,M))
        out[i,5]=action_vertex([off,Kleg,('h',p,None),('h',r,None)],a,m,M)
    return out

def boost(beta):
    g=(1-beta@beta)**-.5
    B=np.zeros((4,4));B[0,0]=g;B[0,1:]=g*beta;B[1:,0]=g*beta
    B[1:,1:]=np.eye(3)+g*g/(g+1)*np.outer(beta,beta)
    return B

def make_gamma():
    sig=[np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1,-1])]
    Z=np.zeros((2,2))
    return np.array([np.diag([1,1,-1,-1]),*[np.block([[Z,s],[-s,Z]]) for s in sig]])


def symbolic():
    z,d0,t0,s,L,u,qe,q2,a,m,M,xi,B=sp.symbols('z d0 t0 s L u qe q2 a m M xi B')
    d=d0+s*L*z
    t=t0+s*u*z
    W=2*a*a*m*m*(2*qe-t-(1-xi)*(q2-M*M)*t/(d-xi*m*m))/(d-m*m)
    chi=2*a*a*M*M*t/(d-xi*m*m)
    pair=sp.factor(W+chi)
    w=1/z
    b=sp.factor(sp.cancel(pair.subs(z,1/z)).subs(z,0))
    assert sp.simplify(b-2*a*a*(M*M-m*m)*u/L)==0
    gp=(xi*m*m-d0)/(s*L)
    R=sp.factor(sp.cancel((z-gp)*pair).subs(z,gp))
    assert sp.simplify(R-2*a*a*q2*t.subs(z,gp)/(s*L))==0
    qQ=(q2+d-M*M)/2
    Vw=4*a*a*m*m*(qe-(1-xi)*qQ*t/(d-xi*m*m))/(d-m*m)
    Vc=-2*a*a*(d-M*M)*t/(d-xi*m*m)
    one=sp.factor(Vw+Vc-W-chi)
    assert sp.simplify(one+2*a*a*t)==0
    t2=qe-t0
    assert sp.simplify(-2*a*a*((t0+u*z)+(t2-u*z))+2*a*a*qe)==0
    return dict(one_ordering_boundary=str(b),unphysical_gauge_pole_residue=str(R),
                one_ordering_ST_difference=str(one),
                total_ST_difference_including_seagull='0',
                total_boundary='2 a m B (q.e) + 2 a^2 (M^2-m^2) (n.e)/(q.n)')

def unit(rng):
    a=rng.normal(size=3)
    return a/np.linalg.norm(a)

def fixture(rng,crossed):
    a=rng.uniform(.4,.75); m=rng.uniform(.75,1.2); M=rng.uniform(.6,1.3)
    B=boost(rng.uniform(.03,.45)*unit(rng))
    ax=unit(rng); x=unit(rng); x-=ax*(x@ax); x/=np.linalg.norm(x); y=np.cross(ax,x)
    k=rng.uniform(.55,1.55); en=np.sqrt(k*k+M*M)
    p=B@np.r_[en,k*ax]; r=B@np.r_[en,-k*ax]; n=B@np.r_[0.,x+1j*y]
    f=rng.uniform(.5,1.8)*np.r_[1.,unit(rng)]; b=rng.uniform(.5,1.8)*np.r_[1.,unit(rng)]
    if crossed:b=-b
    K=f+b
    q=-(K+p+r)
    gam=make_gamma(); sl=lambda v:np.einsum('i,ijk->jk',ETA@v,gam)
    zz=rng.normal(size=4)+1j*rng.normal(size=4); ww=rng.normal(size=4)+1j*rng.normal(size=4)
    e=np.array([zz@sl(b)@g@sl(f)@ww for g in gam]);e/=np.max(np.abs(e))
    return (q,K,p,r,n,e,a,m,M)

def analytic_poles(F,xi):
    q,K,p,r,n,e,a,m,M=F
    L=2*dot(q,n)
    return [(i,s,mu,(mu-dot(q+x,q+x))/(s*L)) for i,(s,x) in enumerate(((1,p),(-1,r))) for mu in sorted(set([m*m,xi*m*m]))]

def contour(center,R,N,fun):
    vals=[]
    for j in range(N):
        z=R*np.exp(2j*np.pi*(j+.317)/N)
        vals.append(fun(center+z,z))
    return np.mean(vals,axis=0)

def cm(x):
    a=np.asarray(x)
    return np.stack((a.real,a.imag),axis=-1).tolist()

def evaluate(seed,N):
    rng=np.random.default_rng(seed);xis=(.43,1.,2.6); maxs={}; out=[]; rejects=0
    def obs(k,v):maxs[k]=max(maxs.get(k,0),float(np.max(np.abs(v))))
    for j in range(N):
        while True:
            F=fixture(rng,bool(j%2))
            q,K,p,r,n,e,a,m,M=F
            allp=[analytic_poles(F,xi) for xi in xis]
            sep=min(abs(x[3]-y[3]) for poles in allp for i,x in enumerate(poles) for y in poles[i+1:])
            near=min(abs(x[3]) for poles in allp for x in poles)
            if abs(dot(q,n))>.18 and near>.16 and sep>.14 and abs(dot(K,K)-m*m)>.18 and abs(dot(q,q)-m*m)>.2:break
            rejects+=1
        H=p+r; v=m/a; Bh=3*M*M/v/(dot(H,H)-M*M)
        expected=2*a*m*Bh*dot(q,e)+2*a*a*(M*M-m*m)*dot(n,e)/dot(q,n)
        row=dict(q=cm(q),K=cm(K),p=cm(p),r=cm(r),n=cm(n),conserved_e=cm(e),a=a,m=m,M=M,expected_boundary=cm(expected),gauges=[])
        obs('normalization_and_shift_constraints',[dot(K,e),dot(n,n),dot(n,p),dot(n,r),dot(p,p)-M*M,dot(r,r)-M*M])
        for xi,poles in zip(xis,allp):
            def fun(z):
                pp=p+z*n; rr=r-z*n
                vals=direct_amplitudes(q,K,pp,rr,e,a,m,M,xi)
                return np.sum(vals,axis=1)*np.array([-1j,m])
            f0=fun(0)
            edge=max(abs(x[3]) for x in poles)
            R=max(2.,3*edge)
            b0=contour(0,R,32,lambda z,_:fun(z))
            b1=contour(0,1.41*R,36,lambda z,_:fun(z))
            scale=max(1,np.max(np.abs(f0)),abs(expected))
            l1=0j; rs=[]
            for k,(_,s,mu,zp) in enumerate(poles):
                radius=.19*min([abs(zp-p2[3]) for i,p2 in enumerate(poles) if i!=k]+[abs(zp),R-abs(zp)])
                res=contour(zp,radius,20,lambda z,w:w*fun(z))
                l1=l1-res/zp
                obs('residue_ST_difference_scaled',(res[0]-res[1])/max(1,np.max(np.abs(res))))
                if xi!=1 and mu==xi*m*m:
                    x=p if s==1 else r
                    Q=q+x+s*zp*n
                    exact=2*a*a*dot(q,q)*dot(Q,e)/(s*2*dot(q,n))
                    obs('unphysical_pole_formula_scaled',(res-exact)/max(1,abs(exact)))
                    obs('nonzero_unphysical_pole_raw',abs(exact))
                rs.append(dict(position=cm(zp),mass_squared=mu,residues=cm(res)))
            obs('finite_ST_scaled',(f0[0]-f0[1])/scale)
            obs('boundary_vs_formula_scaled',(b0-expected)/scale)
            obs('different_outer_contours_scaled',(b0-b1)/scale)
            obs('boundary_necessary_for_reconstruction_scaled',(l1-f0)/scale)
            obs('reconstruction_with_boundary_scaled',(l1+b0-f0)/scale)
            for z in (R*(.31+.44j), R*(-.24+.51j),.17+.29j):
                zz=fun(z)
                obs('shifted_ST_scaled',(zz[0]-zz[1])/max(1,np.max(np.abs(zz))))
                obs('shifted_mass_shell_scaled',[dot(p+z*n,p+z*n)-M*M,dot(r-z*n,r-z*n)-M*M])
            row['gauges'].append(dict(xi=xi,f0=cm(f0),contour_boundary=cm(b0),second_outer_boundary=cm(b1),residues_only_f0=cm(l1),poles=rs))
        out.append(row)
    physical_errors={k:v for k,v in maxs.items() if k not in {'boundary_necessary_for_reconstruction_scaled','nonzero_unphysical_pole_raw'}}
    assert max(physical_errors.values())<5e-8,maxs
    assert min(maxs['boundary_necessary_for_reconstruction_scaled'],maxs['nonzero_unphysical_pole_raw'])>1e-3,maxs
    return dict(accepted=N,rejected=rejects,seed=seed,maxima=maxs,fixtures=out)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--seed',type=int,default=62023);ap.add_argument('--points',type=int,default=8);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    report={"schema":1,"symbolic":symbolic(),"numerical":evaluate(args.seed,args.points),"script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"versions":{"python":platform.python_version(),"numpy":np.__version__,"sympy":sp.__version__}}
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in report['numerical'].items() if k!='fixtures'},indent=2))

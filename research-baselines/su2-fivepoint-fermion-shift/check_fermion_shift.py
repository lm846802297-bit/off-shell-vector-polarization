#!/usr/bin/env python3
"""Fermion-pair Weyl shift, fixed bosonic invariants; see README for scope."""
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

SIG=np.array([np.eye(2),[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],complex)
BARS=SIG*np.array([1.,-1.,-1.,-1.])[:,None,None]
EPS=np.array([[0,1],[-1,0]])
Z0=np.zeros((2,2),complex)
GAMMA=np.array([np.block([[Z0,s],[b,Z0]]) for s,b in zip(SIG,BARS)])
slash=lambda p:np.einsum('i,ijk->jk',ETA@p,GAMMA)
def inv(P):
    return np.array([np.trace(s@P)/2 for s in SIG])
def pvec(l,t):
    return inv(np.outer(l,t))
def cur_outer(l,t):
    return 2*pvec(l,t)

def symbolic():
    z,k,m,M,a,r,s,qe,q2,t,xi=sp.symbols('z k m M a r s A q2 t xi')
    L=sp.Matrix([sp.Symbol('L0'),sp.Symbol('L1')])
    B=sp.Matrix([sp.Symbol('B0'),sp.Symbol('B1')])
    T=sp.Matrix([sp.Symbol('T0'),sp.Symbol('T1')])
    U=sp.Matrix([sp.Symbol('U0'),sp.Symbol('U1')])
    F=(L+z*B)*T.T
    Fb=B*(U-z*T).T
    A=(L+z*B)*(U-z*T).T
    C=B*T.T
    check=lambda mat:all(sp.expand(x)==0 for x in mat)
    assert check(F+Fb-L*T.T-B*U.T)
    assert check(A.diff(z,2)+2*C)
    assert check(A.diff(z,3))
    assert F.det().expand()==Fb.det().expand()==A.det().expand()==C.det().expand()==0
    K=F+Fb
    for J in (A,C):
        # Bilinear K·J is half of tr(adj K J); no complex conjugate in z.
        assert sp.expand(sp.trace(K.adjugate()*J))==0
    d=k-m*m; c=k-xi*m*m
    G=2*a*a*m*m*(2*qe-t-(1-xi)*(q2-M*M)*t/c)/d+2*a*a*M*M*t/c
    diff=sp.factor(G-G.subs(xi,1))
    assert sp.simplify(diff-2*a*a*m*m*(xi-1)*q2*t/(d*c))==0
    return dict(momentum_total_invariant=True,rank_one_null=True,both_currents_conserved=True,
                A_quadratic_coefficient='-j_B',B_holomorphically_constant=True,
                one_ordering_gauge_variation='2 a^2 m^2 (xi-1) q^2 (Q.j)/[(d-m^2)(d-xi m^2)]')

def dirac_currents(lf,tf,lb,tb):
    uR=np.r_[0.,0.,lf]
    uL=np.r_[EPS@tf,0.,0.]
    barR=np.r_[tb,0.,0.]
    barL=np.r_[0.,0.,-lb@EPS]
    return np.array([[r@g@u for g in GAMMA] for r,u in ((barR,uR),(barL,uL))]),(uR,uL),(barR,barL)

def closed_G(F,e,xi):
    q,K,p,r,a,m,M=F
    G=2*a*m*(3*M*M/(m/a))/(dot(p+r,p+r)-M*M)*dot(q,e)
    for x in (p,r):
        Q=q+x; d=dot(Q,Q)
        A=dot(q,e);t=dot(Q,e)
        G+=2*a*a*m*m*(2*A-t-(1-xi)*(dot(q,q)-M*M)*t/(d-xi*m*m))/(d-m*m)+2*a*a*M*M*t/(d-xi*m*m)
    return G

def gauge(F,e,xi):
    q,K,p,r,a,m,M=F
    ans=0j
    for x in (p,r):
        Q=q+x; d=dot(Q,Q)
        ans+=2*a*a*m*m*(xi-1)*dot(q,q)*dot(Q,e)/((d-m*m)*(d-xi*m*m))
    return ans

def grades_from_action(F,spin,z,xi):
    q,K,p,r,a,m,M=F
    lf,tf,lb,tb=spin
    L=lf+z*lb; T=tb-z*tf
    fp,bp=pvec(L,tf),pvec(lb,T)
    J,us,bs=dirac_currents(L,tf,lb,T)
    amp=[]
    for e in J:
        rows=direct_amplitudes(q,K,p,r,e,a,m,M,xi)
        amp.append(np.sum(rows,axis=1)*np.array([-1j,m]))
    return np.array(amp),J,fp,bp,us,bs

def cm(a):
    a=np.asarray(a)
    return np.stack((a.real,a.imag),axis=-1).tolist()

def evaluate(seed,N):
    rng=np.random.default_rng(seed); xis=(.46,1.,2.35); out=[]; rejected=0; errors={}; mins={'fixed_denominator':1e100}
    def obs(k,v): errors[k]=max(errors.get(k,0.),float(np.max(np.abs(v))))
    def spin():return (rng.normal(size=2)+1j*rng.normal(size=2))*.48
    for i in range(N):
        a=rng.uniform(.35,.65);m=rng.uniform(.7,1.3);M=rng.uniform(.6,1.4)
        while True:
            lf,lb=spin(),spin(); tf=lf.conj(); tb=lb.conj()
            if i%2:tb=-tb
            K=pvec(lf,tf)+pvec(lb,tb)
            pv,rv=rng.normal(size=(2,3))
            p=np.r_[np.sqrt(M*M+pv@pv),pv]
            r=np.r_[np.sqrt(M*M+rv@rv),rv]
            q=-(K+p+r); F=q,K,p,r,a,m,M
            gap=min([abs(dot(q+x,q+x)-xi*m*m) for x in (p,r) for xi in xis]+[abs(dot(K,K)-m*m),abs(dot(p+r,p+r)-M*M)])
            if gap>.22 and abs(dot(K,K))>.2 and abs(dot(q,q)-m*m)>.3:break
            rejected+=1
        mins['fixed_denominator']=min(mins['fixed_denominator'],gap)
        JA0=cur_outer(lf,tb)
        JA1=cur_outer(lb,tb)-cur_outer(lf,tf)
        JB=cur_outer(lb,tf)
        coJ=[[JA0,JA1,-JB],[JB,np.zeros(4),np.zeros(4)]]
        row=dict(a=a,m=m,M=M,q=cm(q),K=cm(K),p=cm(p),r=cm(r),external_spinors=cm([lf,tf,lb,tb]),gauges=[])
        for xi in xis:
            co=np.array([[closed_G(F,e,xi) for e in A] for A in coJ])
            res=[];gdiff=[]
            for z in (0j,.31+.77j,-1.15+.48j):
                amp,J,fp,bp,us,bs=grades_from_action(F,(lf,tf,lb,tb),z,xi)
                Jt=np.einsum('nji,j->ni',np.array(coJ),[1,z,z*z])
                expected=co@np.array([1,z,z*z])
                sc=max(1.,np.max(np.abs(amp)),np.max(np.abs(expected)))
                obs('mass_momentum_conservation_and_spinor_equations',[dot(fp,fp),dot(bp,bp),*list(fp+bp-K),*[x for u in us for x in slash(fp)@u],*[x for b in bs for x in b@slash(bp)]])
                obs('current_conservation',J@ETA@K)
                obs('two_component_vs_dirac_current',Jt-J)
                obs('finite_ST_scaled',(amp[:,0]-amp[:,1])/sc)
                obs('direct_vs_closed_scaled',(amp-expected[:,None])/sc)
                if xi!=1:
                    amp1=grades_from_action(F,(lf,tf,lb,tb),z,1.)[0]
                    delta=amp-amp1
                    ee=np.array([gauge(F,e,xi) for e in J])
                    obs('gauge_change_vs_closed_scaled',(delta-ee[:,None])/max(sc,np.max(np.abs(amp1))))
                    obs('nonzero_gauge_change_scaled',delta/max(sc,np.max(np.abs(amp1))))
                    gdiff.append(cm(delta))
                res.append(dict(z=cm(z),amplitudes=cm(amp)))
            R=3.7
            theta=2*np.pi*(np.arange(12)+.283)/12
            zs=R*np.exp(1j*theta)
            vals=np.array([grades_from_action(F,(lf,tf,lb,tb),z,xi)[0] for z in zs])
            four=np.array([np.mean(vals*np.exp(-1j*k*theta)[:,None,None],axis=0)/R**k for k in (-2,-1,0,1,2,3)])
            ex=np.array([co[:,k] for k in range(3)])[:,:,None]
            scale=max(1.,np.max(np.abs(vals)),np.max(np.abs(ex)))
            obs('circle_modes_vs_polynomial_scaled',(four[2:5]-ex)/scale)
            obs('forbidden_fourier_modes_scaled',four[[0,1,5]]/scale)
            obs('flat_current_nonzero_modes_scaled',four[[3,4],1,:]/scale)
            obs('quadratic_minus_flat_link_scaled',(four[4,0,:]+four[2,1,:])/scale)
            obs('nonzero_quadratic_signal_scaled',four[4,0,:]*R**2/scale)
            row['gauges'].append(dict(xi=xi,analytic_coefficients=cm(co),numeric_circle_modes=cm(four),point_checks=res,current_gauge_difference=gdiff))
        out.append(row)
    assert max(v for k,v in errors.items() if not k.startswith('nonzero_'))<5e-8,errors
    assert min(v for k,v in errors.items() if k.startswith('nonzero_'))>1e-3,errors
    return dict(seed=seed,accepted=N,rejected=rejected,minimum=mins,maxima=errors,fixtures=out)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--seed',type=int,default=82023);ap.add_argument('--points',type=int,default=8);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    assert args.points>0
    report={"schema":1,"symbolic":symbolic(),"numerical":evaluate(args.seed,args.points),"script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"versions":{"python":platform.python_version(),"numpy":np.__version__,"sympy":sp.__version__}}
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in report['numerical'].items() if k!='fixtures'},indent=2))

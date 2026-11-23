#!/usr/bin/env python3
"""Test two distinguishable conserved f fbar sources joined by the WW hh tree."""
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

def direct_amplitudes(q,K,p,r,e,a,m,M,xi,root_current):
    """Returns q-contracted vector amplitude and unphased Goldstone amplitude."""
    root=[('W',q,root_current),('c',q,None)]
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
    a,m,d,xi,je,jQ,eQ,jq,M,q2,xi_r=sp.symbols('a m d xi je jQ eQ jq M q2 xi_r')
    D=d-m*m; C=d-xi*m*m
    V=4*a*a*m*m*(je-(1-xi)*jQ*eQ/C)/D
    chi=-2*a*a*(2*jQ-jq)*eQ/C
    ix=sp.factor(V+chi-(V+chi).subs(xi,1))
    assert sp.simplify(ix-2*a*a*m*m*(xi-1)*jq*eQ/(C*D))==0
    sewn=sp.factor((V+chi).subs(jq,0))
    assert sp.simplify(sewn-4*a*a*(m*m*je-jQ*eQ)/D)==0
    root_factor=-(1-xi_r)*jq/(q2-m*m)
    singular=sp.cancel(root_factor.subs(xi_r,q2/(m*m)))
    assert sp.simplify(singular-jq/(m*m))==0
    return dict(internal_sewn_exchange=str(sewn),generic_unconserved_internal_gauge_difference=str(ix),
                root_pole_coefficient_along_q=str(singular),root_pole_if_conserved='0')

def unit(rng):
    x=rng.normal(size=3)
    return x/np.linalg.norm(x)

def fermion_split(P,rng):
    vel=(P[1:]/P[0]).real
    sgn=np.sign(P[0].real)
    k=np.sqrt(dot(P,P).real)
    B=boost(vel); n=unit(rng)
    f,b=[sgn*k/2 * (B@np.r_[1.,s*n]) for s in (1.,-1.)]
    assert np.max(np.abs(f+b-P))<1e-10
    S=make_gamma()
    sl=lambda p:np.einsum('i,ijk->jk',ETA@p,S)
    zz=rng.normal(size=4)+1j*rng.normal(size=4); ww=rng.normal(size=4)+1j*rng.normal(size=4)
    j=np.array([zz@sl(b)@g@sl(f)@ww for g in S])
    j/=np.max(np.abs(j))
    return j,f,b

def build_F(rng):
    a=rng.uniform(.4,.8); m=rng.uniform(.75,1.3); M=rng.uniform(.6,1.15)
    K=rng.uniform(.55,1.8)*(boost(rng.uniform(.03,.5)*unit(rng))@BASIS[0])
    pp,rr=rng.normal(size=(2,3))
    p=np.r_[np.sqrt(M*M+pp@pp),pp]
    r=np.r_[np.sqrt(M*M+rr@rr),rr]
    q=-(K+p+r);z0=dot(q,q)/(m*m)
    eis=[.37,1.,2.8,z0]
    gap=min(abs(dot(q+x,q+x)-xi*m*m) for x in (p,r) for xi in eis)
    if min(gap,abs(dot(K,K)-m*m),abs(dot(q,q)-m*m))<.28:return None
    e,f,b=fermion_split(K,rng)
    j,g,c=fermion_split(q,rng)
    return (q,K,p,r,a,m,M),e,j,(f,b,g,c),float(z0)

def inner(F,E,e,xi):
    q,K,p,r,a,m,M=F
    return -1j*direct_amplitudes(q,K,p,r,e,a,m,M,xi,E)[0]

def oracle(F,j,e):
    q,K,p,r,a,m,M=F; v=m/a; H=p+r
    E=[]
    for x in (p,r):
        Q=q+x
        E.append(4*a*a*(m*m*dot(j,e)-dot(j,Q)*dot(e,Q))/(dot(Q,Q)-m*m))
    return np.array([*E,2*a*m*(3*M*M/v)/(dot(H,H)-M*M)*dot(j,e),2*a*a*dot(j,e)])

def Pprop(p,m,xi):
    return (ETA-(1-xi)*np.outer(p,p)/(dot(p,p)-xi*m*m))/(dot(p,p)-m*m)

def cm(a):
    a=np.asarray(a,dtype=complex);return np.stack((a.real,a.imag),axis=-1).tolist()

def evaluate(seed,N):
    rng=np.random.default_rng(seed);maxs={};cases=[];rejected=0
    def obs(k,v):maxs[k]=max(maxs.get(k,0.),float(np.max(np.abs(v))))
    for _ in range(N):
        F0=build_F(rng)
        while F0 is None:
            rejected+=1;F0=build_F(rng)
        F,e,j,fs,z0=F0
        q,K,p,r,a,m,M=F; Dq=dot(q,q)-m*m;DK=dot(K,K)-m*m
        xi0=(.37,1.,2.8,z0)
        # Explicit action-derived uncontracted WW hh matrices, by graph type.
        # T is bilinear in contravariant arguments E,e.
        matrices={xi:np.array([[inner(F,E,u,xi) for u in BASIS] for E in BASIS]).transpose(2,0,1) for xi in xi0}
        first=matrices[1.]
        e2=oracle(F,j,e)
        expected=sum(e2)/(Dq*DK)
        scale=max(1.,abs(expected),np.max(np.abs(e2)))
        obs('physical_sources_on_shell_and_transverse',[dot(f,f) for f in fs]+[dot(q,j),dot(K,e)])
        fcut=lambda X,j0,e0:np.array([j0@A@e0 for A in X])
        def attached(j0,xi_root,xi_far,xi_inner, subtract_roundoff=False):
            qh=np.asarray(q,dtype=np.longdouble)
            kh=np.asarray(K,dtype=np.longdouble)
            mh=np.longdouble(m)
            xh=np.longdouble(xi_root)
            Jh=np.asarray(j0,dtype=np.clongdouble)
            Jr=Jh@ETA@Pprop(qh,mh,xh)  # contravariant after propagation
            if subtract_roundoff:
                # Test the *same* tensor with its explicitly measured
                # longitudinal floating-point source contamination removed.
                dqh=dot(qh,qh)-mh*mh; cq=dot(qh,qh)-xh*mh*mh
                Jr=Jr+(1-xh)*dot(Jh,qh)*qh/(dqh*cq)
            er=np.asarray(e,dtype=np.clongdouble)@ETA@Pprop(kh,mh,np.longdouble(xi_far))
            return sum(fcut(matrices[xi_inner],Jr,er))
        rec=dict(q=cm(q),K=cm(K),p=cm(p),r=cm(r),gauge_root_value=z0,a=a,m=m,M=M,
                 source_near=cm(j),source_far=cm(e),external_fermions=cm(fs),closed_graph_groups=cm(e2),regular=[],near_root=[])
        for xi in xi0:
            val=fcut(matrices[xi],j,e)
            regroup=np.r_[val[:2].sum(),val[2:4].sum(),val[4:]]
            obs('internal_exchange_vs_closed_scaled',(regroup-e2)/scale)
            if xi!=1:
                obs('missing_chi_negative_control_scaled',(val-fcut(first,j,e))[[0,2]]/scale)
                wrong=j+.19*q
                exact=[]
                for x in (p,r):
                    Q=q+x;D=dot(Q,Q)-m*m;C=dot(Q,Q)-xi*m*m
                    exact.append(2*a*a*m*m*(xi-1)*dot(q,wrong)*dot(e,Q)/(C*D))
                diff=fcut(matrices[xi]-first,wrong,e)
                targ=np.array([diff[:2].sum(),diff[2:4].sum()])
                obs('unconserved_internal_gauge_formula_scaled',(targ-exact)/max(scale,np.max(np.abs(targ))))
                obs('unconserved_internal_negative_control_scaled',targ/scale)
            if xi in (.37,1.,2.8):
                v=attached(j,xi,.71,xi)
                w=attached(j,.71,xi,xi)
                obs('independent_outer_gauges_scaled',[abs(v-expected)/scale,abs(w-expected)/scale])
                rec['regular'].append(dict(xi=xi,by_graph=cm(val),near_gauge_varied=cm(v),far_gauge_varied=cm(w)))
        # Vary the near propagator through c=q^2-xi_root*m^2=0 at *fixed*
        # physical fermion momenta. The interior is analytic at xi_inner=z0.
        T=matrices[z0]
        wrong=j+.19*q
        pole=(dot(q,wrong)/(m*m))*sum(fcut(T,q,e))/DK
        for gap in (1e-3,1e-5,1e-7):
            xm,xp=z0-gap/(m*m),z0+gap/(m*m)
            vp,vm=attached(j,xm,.71,z0),attached(j,xp,.71,z0)
            fp,fm=attached(wrong,xm,.71,z0),attached(wrong,xp,.71,z0)
            res=gap*(fp-fm)/2
            cp,cmn=attached(j,xm,.71,z0,True),attached(j,xp,.71,z0,True)
            obs('root_raw_roundoff_diagnostic',[abs(vp-expected)/scale,abs(vm-expected)/scale])
            obs('conserved_root_pole_probe_scaled',[abs(cp-expected)/scale,abs(cmn-expected)/scale])
            obs('unconserved_root_residue_vs_exact_scaled',(res-pole)/max(1,np.abs(pole)))
            obs('unconserved_root_pole_negative_control_scaled',res/scale)
            rec['near_root'].append(dict(gap=gap,conserved_raw=cm([vp,vm]),conserved_corrected=cm([cp,cmn]),nonconserved=cm([fp,fm]),residue=cm(res),expected_wrong_residue=cm(pole)))
        # Root pole numerator projected with a conserved current; unitary/physical numerator.
        obs('sewn_physical_pole_numerator',j@ETA@((ETA-np.outer(q,q)/(m*m))-ETA))
        cases.append(rec)
    assert max(v for k,v in maxs.items() if 'negative_control' not in k and k!='root_raw_roundoff_diagnostic')<5e-7,maxs
    assert min(v for k,v in maxs.items() if 'negative_control' in k)>1e-3,maxs
    return dict(accepted=N,rejected=rejected,seed=seed,maxima=maxs,fixtures=cases)

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--seed',type=int,default=102023); ap.add_argument('--points',type=int,default=8); ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    assert args.points>0
    data={"schema":1,"symbolic":symbolic(),"numerical":evaluate(args.seed,args.points),"script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"versions":{"python":platform.python_version(),"numpy":np.__version__,"sympy":sp.__version__}}
    args.output.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in data['numerical'].items() if k!='fixtures'},indent=2))

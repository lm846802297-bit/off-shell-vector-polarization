#!/usr/bin/env python3
"""Same-line double Higgs and virtual h splitting: signed neutral-sector tree checks."""
import argparse
import hashlib
import itertools
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

def boost(beta):
    g=(1-beta@beta)**-.5
    B=np.zeros((4,4));B[0,0]=g;B[0,1:]=g*beta;B[1:,0]=g*beta
    B[1:,1:]=np.eye(3)+g*g/(g+1)*np.outer(beta,beta)
    return B

def make_gamma():
    sig=[np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1,-1])]
    Z=np.zeros((2,2))
    return np.array([np.diag([1,1,-1,-1]),*[np.block([[Z,s],[-s,Z]]) for s in sig]])




GAM=make_gamma()
G5=1j*GAM[0]@GAM[1]@GAM[2]@GAM[3]
I4=np.eye(4)
PL=(I4-G5)/2
def slash(k):return np.einsum('i,ijk->jk',ETA@k,GAM)
def unit(rng):
    x=rng.normal(size=3)
    return x/np.linalg.norm(x)
def cm(x):
    a=np.asarray(x,dtype=complex)
    return np.stack((a.real,a.imag),axis=-1).tolist()

def get_fixture(rng):
    a=rng.uniform(.4,.8);m=rng.uniform(.75,1.25);M=rng.uniform(.6,1.35);mf=rng.uniform(.18,.69)
    K=-rng.uniform(.6,1.75)*(boost(rng.uniform(.05,.5)*unit(rng))@BASIS[0])
    d=rng.normal(size=(2,3))
    p,r=[np.r_[-np.sqrt(M*M+v@v),v] for v in d]
    q=-(K+p+r);s=dot(q,q);D=s-m*m;xr=float(s/m**2)
    if min([abs(D),abs(dot(K,K)-m*m)]+[abs(dot(q+x,q+x)-xi*m*m) for x in (p,r) for xi in (.38,1.,2.5,xr)])<.3:return
    if s<=4*mf*mf:return
    def spin_pair(P,mf):
        E=np.sqrt(dot(P,P))/2
        pp=np.sqrt(E**2-mf**2)*unit(rng)
        L=boost(P[1:]/P[0]); sg=np.sign(P[0])
        f,b=[sg*L@np.r_[E,k*pp] for k in (1.,-1.)]
        rr=rng.normal(size=4)+1j*rng.normal(size=4)
        ll=rng.normal(size=4)+1j*rng.normal(size=4)
        us=(slash(f)+mf*I4)@rr; bs=ll@(slash(b)-mf*I4)
        us/=np.linalg.norm(us); bs/=np.linalg.norm(bs)
        return f,b,us,bs
    f,b,u,vb=spin_pair(q,mf)
    fk,bk,uk,vbk=spin_pair(K,0.)
    e=np.array([vbk@g@uk for g in GAM])
    return (a,m,M,mf,q,K,p,r,xr),e,(f,b,u,vb),(fk,bk,uk,vbk)

def source_vertices(a,m,mf):
    # From i times the field coefficient of the kinetic term for T3=-1/2,
    # and from -y[(v+h+i chi)bar f_L f_R + (v+h-i chi)bar f_R f_L]/sqrt(2).
    v=m/a; PR=I4-PL
    Lyuk=lambda h,chi:-(mf/v)*((v+h+1j*chi)*PR+(v+h-1j*chi)*PL)
    chi=1j*(Lyuk(0,1)-Lyuk(0,-1))/2
    return np.array([-1j*a*g@PL for g in GAM]),chi



def symbolic():
    x,y,z,L,R,PLs,PRs,G5s,mu=sp.symbols('x y z L R PL PR gamma5 mu',commutative=False)
    contacts=y*PLs+PLs*x-y*PRs-PRs*x
    assert sp.expand(contacts.subs(PRs,PLs+G5s)+y*G5s+G5s*x)==0
    d,C,xi,m,U,J=sp.symbols('d C xi m U J')
    sewn=J/d-(1-xi)*U/(d*C)
    assert sp.simplify(sewn.subs(U,0)-J/d)==0
    assert sp.simplify(sewn-sewn.subs(xi,1)-(xi-1)*U/(d*C))==0
    return dict(fixed_Higgs_order_endpoint='-S(f+x+K) gamma5 - gamma5 S(f+x)',
                splitting_endpoint='-gamma5',outer_gauge_delta='(xi-1)(K.J)(K.e)/(D_K C_K)')

def ordered_chain(order,vertex,fs,u,vbar,SP):
    # Independent raw route: compose bare matrices in oriented absorption order.
    L=fs.copy(); st=u.copy()
    for i,(s,v) in enumerate(order):
        st=vertex(s)@st
        L=L+v
        if i+1<len(order): st=(1j*SP(L))@st
    return vbar@st

def evaluate(seed,N):
    rng=np.random.default_rng(seed)
    xis=(.43,1.,2.65); cases=[]; rejected=0; maxs={}
    def obs(k,x):maxs[k]=max(maxs.get(k,0.),float(np.max(np.abs(x))))
    for _ in range(N):
        while True:
            f0=get_fixture(rng)
            if f0 is not None:
                (a,m,M,mu,q,K,p,r,xr),e,(f,g,u,vb),far=f0
                H=p+r
                ks=(p,r,K)
                paths=[f+x for x in ks]+[f+x+y for x,y in itertools.permutations(ks,2)]+[f+H]
                gap=min([abs(dot(l,l)-mu*mu) for l in paths]+[abs(dot(K,K)-xi*m*m) for xi in xis]+[abs(dot(H,H)-M*M)])
                if gap>.22: break
            rejected+=1
        v=m/a;b=mu/v;D=dot(K,K)-m*m; DH=dot(H,H)-M*M
        W,_=source_vertices(a,m,mu)
        Y=lambda h: -mu*(1+h/v)*I4
        Hvertex=1j*(Y(1)-Y(-1))/2
        hhh=action_vertex([('h',p,None),('h',r,None),('h',-H,None)],a,m,M)
        gamma=GAM@PL
        SP=lambda l:(slash(l)+mu*I4)/(dot(l,l)-mu*mu)
        P=vb@G5@u
        pref=3*M*M/(v*DH)
        rec=dict(a=a,m=m,M=M,fermion_mass=mu,q=cm(q),K=cm(K),p=cm(p),r=cm(r),massive_pair=cm([f,g]),far_massless_pair=cm(far[:2]),minimum_denominator=gap,gauges=[])
        obs('on_shell_and_fermion_equations',[dot(f,f)-mu*mu,dot(g,g)-mu*mu,dot(K,e),*list((slash(f)-mu*I4)@u),*list(vb@(slash(g)+mu*I4))])
        obs('action_scalar_splitting_rule',hhh+1j*3*M*M/v)
        near_raw=[]; near_expected=[]
        for x,y in ((p,r),(r,p)):
            # Materialize the grouped closed form from explicit momenta.
            grouped=lambda B: np.array([[vb@SP(f+K+x)@SP(f+K)@h@u, vb@SP(f+x+K)@h@SP(f+x)@u, vb@h@SP(f+x+y)@SP(f+x)@u] for h in B]).T
            j=grouped(gamma)
            c=grouped([G5]).sum()
            rem=vb@(SP(f+x+K)@G5+G5@SP(f+x))@u
            obs('ordered_double_H_endpoint_ward',(j.sum(axis=0)@ETA@K+mu*c+rem)/max(1.,abs(rem),abs(mu*c)))
            parts=[]
            for order in [(('W',K),('h',x),('h',y)),(('h',x),('W',K),('h',y)),(('h',x),('h',y),('W',K))]:
                parts.append([ordered_chain(order,lambda typ:ww if typ=='W' else Hvertex,f,u,vb,SP) for ww in W])
            near_raw.extend(parts); near_expected.extend(b*b*j)
        jH=np.array([[vb@SP(f+K)@h@u,vb@h@SP(f+H)@u] for h in gamma]).T
        pH=vb@(SP(f+K)@G5+G5@SP(f+H))@u
        obs('virtual_H_endpoint_ward',(jH.sum(axis=0)@ETA@K+mu*pH+P)/max(1.,abs(P),abs(mu*pH)))
        for order in [(('W',K),('h',H)),(('h',H),('W',K))]:
            near_raw.append([ordered_chain(order,lambda typ:ww if typ=='W' else Hvertex,f,u,vb,SP)*(1j/DH)*hhh for ww in W])
        near_expected.extend(b*pref*jH)
        near_raw=np.array(near_raw); near_expected=np.array(near_expected)
        assert near_raw.shape==near_expected.shape==(8,4)
        obs('raw_vs_closed_near_sources_scaled',(near_raw/(-1j*a)-near_expected)/max(1.,np.max(np.abs(near_expected))))
        def two(vs):
            return np.array([vs[:6].sum(),vs[6:].sum()])
        def sew(E,xi):
            return near_raw@ETA@vector_prop(K,m,xi)@ETA@E/(-a)
        pred=near_expected@ETA@e/D
        true_b=two(sew(e,1.)); ec=e+.17*K;bad_b=two(sew(ec,1.))
        jclasses=np.array([near_expected[:6].sum(axis=0),near_expected[6:].sum(axis=0)])
        phys_scale=np.maximum(1.,np.abs(two(pred)))
        obs('nonzero_physical_double_H_amplitude',true_b[0]/phys_scale[0])
        obs('nonzero_physical_splitting_amplitude',true_b[1]/phys_scale[1])
        # Independently re-enumerate permutations with swapped label dictionary.
        pa,pb=(r,p); orders=[]
        for x,y in ((pa,pb),(pb,pa)):
            for j in range(3):
                order=[('h',x),('h',y)];order.insert(j,('W',K));orders.append(order)
        bose=np.array([[ordered_chain(order,lambda typ:ww if typ=='W' else Hvertex,f,u,vb,SP) for ww in W] for order in orders])
        obs('Bose_order_block_exchange_scaled',(bose-near_raw[np.r_[3:6,0:3]])/max(1.,np.max(np.abs(near_raw))))
        for xi in xis:
            got=sew(e,xi); bad=two(sew(ec,xi))
            C=dot(K,K)-xi*m*m
            target=(xi-1)*(jclasses@ETA@K)*dot(K,ec)/(D*C)
            obs('physical_gauge_variation_scaled',(two(got)-true_b)/phys_scale)
            obs('each_order_direct_vs_closed_scaled',(got-pred)/max(1.,np.max(np.abs(pred))))
            diff=bad-bad_b
            bscale=np.maximum(1.,np.maximum(np.abs(bad),np.abs(bad_b)))
            obs('wrong_source_analytic_variation_scaled',(diff-target)/bscale)
            obs('nonzero_wrong_source_double_H',diff[0]/bscale[0])
            obs('nonzero_wrong_source_splitting',diff[1]/bscale[1])
            rec['gauges'].append(dict(xi=xi,raw_six_plus_two=cm(got),closed_six_plus_two=cm(pred),physical_by_class=cm(two(got)),
                                      wrong_source_var_by_class=cm(diff),wrong_source_expectation=cm(target)))
        cases.append(rec)
    assert max(v for k,v in maxs.items() if not k.startswith('nonzero_'))<1e-9,maxs
    assert min(v for k,v in maxs.items() if k.startswith('nonzero_'))>1e-4,maxs
    return dict(accepted=N,rejected=rejected,seed=seed,maxima=maxs,fixtures=cases)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--seed',type=int,default=132023);ap.add_argument('--points',type=int,default=8);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    assert args.points>0
    data={'schema':1,'symbolic':symbolic(),'numerical':evaluate(args.seed,args.points),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'versions':{'python':platform.python_version(),'numpy':np.__version__,'sympy':sp.__version__}}
    args.output.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in data['numerical'].items() if k!='fixtures'},indent=2))

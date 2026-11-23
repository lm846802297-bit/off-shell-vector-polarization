#!/usr/bin/env python3
"""Physical chiral Yukawa source sewn to the earlier SU(2) WW hh subtree."""
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
    u,T,V,G,m,D,C,xi,c=sp.symbols('u T V G m D C xi c')
    root=(T-(1-xi)*u*V/C)/D-u*G/(m*C)
    target=(T-u*V/m**2)/D
    assert sp.simplify((root-target).subs({G:V/m,C:D+(1-xi)*m*m}))==0
    delta=(c*u-c*u*(D+m*m)/m**2)/D
    assert sp.simplify(delta+c*u/m**2)==0
    return dict(convention_source_current='q.j=m_f vbar gamma5 u_f',
                root_cancellation='(j.T - (q.j) V/m^2)/(q^2-m^2)',
                inner_gauge_variation_after_root_sewing='-c(q.j)/m^2',
                c_per_ordering='2 a^2 m^2 (xi-1) (e.Q) / ((Q^2-m^2)(Q^2-xi m^2))')

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

def evaluate(seed,N):
    rng=np.random.default_rng(seed);maxs={};cases=[];rejected=0
    def obs(k,v):maxs[k]=max(maxs.get(k,0.),float(np.max(np.abs(v))))
    for _ in range(N):
        F=get_fixture(rng)
        while F is None:
            rejected+=1;F=get_fixture(rng)
        (a,m,M,mf,q,K,p,r,xr),e,(f,b,u,vb),(fk,bk,uk,vbk)=F
        D=dot(q,q)-m*m
        j=np.array([vb@g@PL@u for g in GAM]); Wv,chv=source_vertices(a,m,mf)
        Jraw=np.array([vb@g@u for g in Wv]); Craw=vb@chv@u
        jq=dot(q,j); Ps=vb@G5@u
        obs('spinor_mass_shell_and_current_conservation',[dot(f,f)-mf*mf,dot(b,b)-mf*mf,dot(fk,fk),dot(bk,bk),dot(K,e),*list(f+b-q),*list(fk+bk-K)])
        obs('spinor_dirac_equations',[*list((slash(f)-mf*I4)@u),*list(vb@(slash(b)+mf*I4))])
        obs('source_signed_Ward_residual',[jq-mf*Ps,Craw-a*jq/m,*list(Jraw+1j*a*j)])
        obs('nonzero_massive_source_divergence',jq)
        # Independent tensor/current route: coefficient extraction at the right
        # vertices, no Ward projection or replacement in the raw sewing.
        xi_set=(.38,1.,2.5,xr)
        T={}
        for xi in xi_set:
            raw=np.array([direct_amplitudes(q,K,p,r,e,a,m,M,xi,z) for z in BASIS])
            T[xi]=(-1j*np.sum(raw[:,0,:],axis=1),np.sum(raw[0,1,:]))
        def sew(t,chi,x):
            C=dot(q,q)-x*m*m
            Prop=(ETA-(1-x)*np.outer(q,q)/C)/D
            ww=(Jraw@ETA@(-1j*Prop)@(1j*t))/(-1j*a)
            cc=(Craw*(1j/C)*chi)/(-1j*a)
            return ww,cc
        t1,g1=T[1.]
        J_eff=j-jq*q/(m*m)
        oracle_1=J_eff@t1/D
        scale=max(1.,abs(oracle_1))
        rec=dict(m=m,M=M,mf=mf,a=a,q=cm(q),K=cm(K),p=cm(p),r=cm(r),massive_pair=cm([f,b]),
                 massless_pair=cm([fk,bk]),current_massive=cm(j),current_massless=cm(e),source_divergence=cm(jq),source_chi_raw=cm(Craw),xi_at_root=xr,regular=[],pole_probes=[])
        for xi in xi_set:
            t,g=T[xi]
            obs('right_ST_relation_scaled',(q@t-m*g)/max(1.,np.max(np.abs(t)),abs(m*g)))
            cs=0j
            for x in (p,r):
                Q=q+x;d=dot(Q,Q)
                cs+=2*a*a*m*m*(xi-1)*dot(e,Q)/((d-m*m)*(d-xi*m*m))
            expected=-jq*cs/(m*m)
            aa=J_eff@t/D
            obs('inner_residual_formula_scaled',(aa-oracle_1-expected)/max(scale,abs(aa),abs(expected)))
            obs('nonzero_inner_residual_scaled',(aa-oracle_1)/max(scale,abs(aa)))
            # Vary root gauge at fixed inner gauge, then compare same-gauge
            # cases for all choices for which the roots are away from zero.
            br={}
            for root_xi in (.61,1.3,2.8):
                vec,gold=sew(t,g,root_xi)
                obs('root_gauges_scaled',(vec+gold-aa)/max(scale,abs(aa)))
                br[str(root_xi)]=cm([vec,gold])
            if abs(dot(q,q)-xi*m*m)>.3:
                vw,cw=sew(t,g,xi)
                obs('same_gauge_root_pair_scaled',(vw+cw-aa)/max(scale,abs(aa)))
                br['same_xi']=cm([vw,cw])
            rec['regular'].append(dict(xi_inner=xi,projection=cm(aa),expected_inner_delta=cm(expected),direct_root_pieces=br))
        tx,gx=T[xr];aa=J_eff@tx/D
        pole=jq*(q@tx)/(m*m)
        obs('nonzero_root_Wonly_residue_scaled',pole/max(scale,abs(aa)))
        for gap in (1e-3,1e-5,1e-6):
            wp,gp=sew(tx,gx,xr-gap/(m*m));wn,gn=sew(tx,gx,xr+gap/(m*m))
            Rw=gap*(wp-wn)/2;Rc=gap*(gp-gn)/2
            errscale=max(scale,abs(aa))
            obs('near_root_full_scaled',[abs(wp+gp-aa)/errscale,abs(wn+gn-aa)/errscale])
            obs('root_split_residue_formula_scaled',[abs(Rw-pole)/max(1.,abs(pole)),abs(Rc+pole)/max(1.,abs(pole))])
            obs('root_total_residue_scaled',(Rw+Rc)/max(1.,abs(pole)))
            rec['pole_probes'].append(dict(gap=gap,W_pair=cm([wp,wn]),chi_pair=cm([gp,gn]),W_residue=cm(Rw),chi_residue=cm(Rc),predicted_W_residue=cm(pole)))
        cases.append(rec)
    assert max(v for k,v in maxs.items() if not k.startswith('nonzero_'))<1e-6,maxs
    assert min(v for k,v in maxs.items() if k.startswith('nonzero_'))>1e-4,maxs
    return dict(accepted=N,rejected=rejected,seed=seed,maxima=maxs,fixtures=cases)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--seed',type=int,default=112023);ap.add_argument('--points',type=int,default=8);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    assert args.points>0
    data={'schema':1,'symbolic':symbolic(),'numerical':evaluate(args.seed,args.points),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'versions':{'python':platform.python_version(),'numpy':np.__version__,'sympy':sp.__version__}}
    args.output.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in data['numerical'].items() if k!='fixtures'},indent=2))

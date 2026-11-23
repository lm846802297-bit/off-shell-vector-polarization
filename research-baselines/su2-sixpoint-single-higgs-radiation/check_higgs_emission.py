#!/usr/bin/env python3
"""One Higgs from the massive line, one from the W/chi bridge."""
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


def current_identity():
    P,Ph,mu,V,G,m,b,D,C,xi,u,eQ,a=sp.symbols('P Ph mu V G m b D C xi u eQ a')
    L=mu*Ph+P
    raw=-b*((1-xi)*L*V/(D*C)+mu*Ph*G/(m*C))
    rhs=-b*L*V/(m*m*D)+b*P*V/(m*m*C)
    assert sp.simplify((raw-rhs).subs({G:V/m,C:D+(1-xi)*m*m}))==0
    variation=b*P*V/(m*m)*(1/C-1/D)
    expected=2*a*a*(xi-1)*u*eQ/(D*C)
    rule={V:2*a*m*eQ,b:a*mu/m,P:u/mu}
    assert sp.simplify((variation.subs(rule)-expected).subs(C,D+(1-xi)*m*m))==0
    return dict(ordered_fermion_endpoint_remainder='Q.J_h-mu P_h = vbar gamma5 u',
                radiative_gauge_change='2 a^2 (xi-1) (q.j) (Q.e)/(D_Q C_Q)',
                paired_exchange_plus_yukawa_gauge_variation='0')

def evaluate(seed,N):
    rng=np.random.default_rng(seed);maxs={};cases=[];rejected=0
    def obs(k,v):maxs[k]=max(maxs.get(k,0.),float(np.max(np.abs(v))))
    gauges=(.43,1.,2.65)
    for _ in range(N):
        while True:
            F=get_fixture(rng)
            if F is not None:
                (a,m,M,mu,q,K,p,r,xr),e,(f,g,u,vb),far=F
                dmin=min([abs(dot(t,t)-mu*mu) for x in (p,r) for t in (f+x,-g-x)]+[abs(dot(q+x,q+x)-xi*m*m) for x in (p,r) for xi in gauges])
                if dmin>.22:break
            rejected+=1
        h=mu/(m/a);W,Chi=source_vertices(a,m,mu)
        # Direct derivative of the real Yukawa H field; no Ward substitution.
        MC=lambda hval: -mu*(1+hval/(m/a))*I4
        H=1j*(MC(1)-MC(-1))/2
        P=vb@G5@u;j=np.array([vb@gg@PL@u for gg in GAM]);dq=dot(q,j)
        obs('on_shell_and_source_ward',[dot(f,f)-mu*mu,dot(g,g)-mu*mu,dot(K,e),dq-mu*P])
        obs('external_dirac_equations',[(slash(f)-mu*I4)@u,vb@(slash(g)+mu*I4)])
        row=dict(a=a,m=m,M=M,mu=mu,q=cm(q),K=cm(K),p=cm(p),r=cm(r),massive_pair=cm([f,g]),
                 far_massless_pair=cm(far[:2]),source_divergence=cm(dq),source_pseudoscalar=cm(P),assignments=[])
        for i,(x,y) in enumerate(((p,r),(r,p))):
            Q=q+x;F1=f+x;F2=-g-x
            SP=lambda t:(slash(t)+mu*I4)/(dot(t,t)-mu*mu)
            S1=SP(F1);S2=SP(F2)
            JL=np.array([vb@gg@PL@S1@u for gg in GAM]);JR=np.array([vb@S2@gg@PL@u for gg in GAM])
            CL=vb@G5@S1@u;CR=vb@S2@G5@u
            J=JL+JR;C=CL+CR;D=dot(Q,Q)-m*m
            Lw=np.array([vb@v@(1j*S1)@H@u for v in W])
            Rw=np.array([vb@H@(1j*S2)@v@u for v in W])
            Lc=vb@Chi@(1j*S1)@H@u
            Rc=vb@H@(1j*S2)@Chi@u
            Vfar=np.array([action_vertex([('W',Q,z),('W',K,e),('h',y,None)],a,m,M) for z in BASIS])
            Cfar=action_vertex([('c',Q,None),('W',K,e),('h',y,None)],a,m,M)
            T=-1j*Vfar; V=Q@T
            obs('far_vertex_ward_and_action',[V-m*Cfar,*list(T-(2*a*m)*(ETA@e))])
            obs('both_ordered_near_source_operators',[np.max(np.abs(Lw+1j*a*h*JL)),np.max(np.abs(Rw+1j*a*h*JR)),abs(Lc-h*h*CL),abs(Rc-h*h*CR)])
            obs('radiative_source_ward_identity_scaled',(dot(Q,J)-mu*C-P)/max(1.,abs(P),abs(dot(Q,J)),abs(mu*C)))
            obs('endpoint_separate_residues_scaled',[dot(Q,JL)-mu*CL-vb@(I4-PL)@u,dot(Q,JR)-mu*CR+vb@PL@u])
            def raw(xi):
                prop=vector_prop(Q,m,xi)
                chprop=1j/(dot(Q,Q)-xi*m*m)
                ww=lambda z:(z@ETA@prop@Vfar)/(-1j*a)
                cc=lambda z:(z*chprop*Cfar)/(-1j*a)
                return np.array([[ww(Lw),cc(Lc)],[ww(Rw),cc(Rc)]])
            def oldraw(xi):
                # Companion shift of the *one targeted exchange ordering only*:
                # prior root-sewn expression, with its physical Goldstone partner.
                raw_lines=np.array([direct_amplitudes(q,K,p,r,e,a,m,M,xi,z) for z in BASIS])
                t=-1j*np.sum(raw_lines[:,0,2*i:2*i+2],axis=1)
                return (j-dq*q/(m*m))@t/(dot(q,q)-m*m)
            base=raw(1.);old1=oldraw(1.)
            asg=dict(label=i,Higgs_on_fermion=cm(x),Higgs_on_bridge=cm(y),Q=cm(Q),near_internal_momenta=cm([F1,F2]),
                     J_h=cm(J),P_h=cm(C),gauges=[])
            for xi in gauges:
                Cq=dot(Q,Q)-xi*m*m
                expected=h*(J@T-dot(Q,J)*V/(m*m))/D+h*P*V/(m*m*Cq)
                parts=raw(xi);Fnew=parts.sum();Fprev=oldraw(xi)
                delta=2*a*a*(xi-1)*dq*dot(e,Q)/(D*Cq)
                scale=max(1.,abs(expected),abs(Fnew),abs(Fprev),abs(old1),np.max(np.abs(parts)),np.max(np.abs(base)))
                dn=Fnew-base.sum();dp=Fprev-old1
                obs('direct_product_vs_factorized_scaled',(Fnew-expected)/scale)
                obs('new_gauge_variation_vs_exact_scaled',(dn-delta)/scale)
                obs('prior_targeted_gauge_variation_vs_exact_scaled',(dp+delta)/scale)
                obs('combined_paired_gauge_variation_scaled',(dn+dp)/scale)
                obs('nonzero_separate_variation_scaled',np.array([dn,dp])/scale)
                obs('nonzero_omitting_one_fermion_order_scaled',((parts-base)[0,:].sum()+dp)/scale)
                obs('nonzero_omitting_Goldstone_bridge_scaled',((parts-base)[:,0].sum()+dp)/scale)
                asg['gauges'].append(dict(xi=xi,new_ordered_near_source_W_chi=cm(parts),old_paired_exchange=cm(Fprev),
                                          new_delta_xi=cm(dn),old_delta_xi=cm(dp),analytic_new_delta=cm(delta)))
            row['assignments'].append(asg)
        cases.append(row)
    assert max(v for k,v in maxs.items() if not k.startswith('nonzero_'))<1e-9,maxs
    assert min(v for k,v in maxs.items() if k.startswith('nonzero_'))>1e-4,maxs
    return dict(accepted=N,rejected=rejected,seed=seed,maxima=maxs,fixtures=cases)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--seed',type=int,default=122023);ap.add_argument('--points',type=int,default=8);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    assert args.points>0
    data={'schema':1,'symbolic':current_identity(),'numerical':evaluate(args.seed,args.points),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'versions':{'python':platform.python_version(),'numpy':np.__version__,'sympy':sp.__version__}}
    args.output.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in data['numerical'].items() if k!='fixtures'},indent=2))

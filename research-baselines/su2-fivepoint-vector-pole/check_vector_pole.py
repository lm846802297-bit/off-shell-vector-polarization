#!/usr/bin/env python3
"""One ff|Whh W-pole check: algebra plus direct Fourier-action vertices."""
import argparse
import hashlib
import json
import platform
from pathlib import Path
import numpy as np
import sympy as sp

ETA = np.diag([1., -1., -1., -1.])
BASIS = np.eye(4)

def symbolic():
    a,m,M,xi,d,Qe,q_e,p_e,r_e,qQ,pQ = sp.symbols('a m M xi Q2 Qe qe pe re qQ pQ')
    c=(1-xi)/(d-xi*m*m)
    # All entries below are reduced, i.e. the raw W-root contractions are i times qv.
    qvW=4*a*a*m*m*(q_e-c*qQ*Qe)/(d-m*m)
    mG_W=2*a*a*m*m*(q_e-p_e-c*(qQ-pQ)*Qe)/(d-m*m)
    qvChi=-a*a*(d-M*M)*2*Qe/(d-xi*m*m)
    mG_Chi=a*a*M*M*2*Qe/(d-xi*m*m)
    vdiff=sp.factor((qvW-mG_W).subs({p_e:Qe-q_e,pQ:d-qQ}))
    gdiff=sp.factor(qvChi-mG_Chi)
    pair = sp.factor(vdiff+gdiff)
    assert sp.expand(pair+2*a*a*Qe) == 0
    extra = sp.expand(pair.subs(Qe,q_e+p_e)+pair.subs(Qe,q_e+r_e)+2*a*a*q_e)
    assert sp.simplify(extra.subs(r_e,-q_e-p_e)) == 0
    # Scalar-exchange term: h momentum H=p+r=-(q+K), with K·e=0.
    H_e, B = sp.symbols('H_e B')
    hs = a*m*B * (2*q_e-(q_e-H_e))
    assert sp.simplify(hs.subs(H_e,-q_e)) == 0
    assert sp.simplify(vdiff+gdiff).has(xi) is False
    return dict(W_exchange=str(vdiff),chi_exchange=str(gdiff),
                one_pair=str(pair),two_pairs_plus_contact="0",scalar_exchange="0")


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
        W=np.zeros(4); h=0.; chi=0.; dh=np.zeros(4,dtype=complex); dc=dh.copy()
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

def analytic_amplitudes(q,K,p,r,e,a,m,M,xi):
    qv=[]; gc=[]; v=m/a
    for x,y in ((p,r),(r,p)):
        Q=q+x; Qe=dot(Q,e); d=dot(Q,Q); f=(1-xi)/(d-xi*m*m)
        qv+= [(2*a*m)**2*(dot(q,e)-f*dot(q,Q)*Qe)/(d-m*m),
               -a*a*dot(q,Q+x)*dot(Q-y,e)/(d-xi*m*m)]
        gc+= [2*a*a*m*(dot(q-x,e)-f*dot(q-x,Q)*Qe)/(d-m*m),
                a*M*M/v*dot(Q-y,e)/(d-xi*m*m)]
    H=p+r; B=3*M*M/v/(dot(H,H)-M*M)
    return np.array([1j*np.array(qv+[2*a*m*B*dot(q,e),2*a*a*dot(q,e)]),
                      gc+[a*B*dot(q-H,e),0]])

def boost(beta):
    g=(1-beta@beta)**-.5
    B=np.zeros((4,4));B[0,0]=g;B[0,1:]=g*beta;B[1:,0]=g*beta
    B[1:,1:]=np.eye(3)+g*g/(g+1)*np.outer(beta,beta)
    return B

def make_gamma():
    sig=[np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1,-1])]
    Z=np.zeros((2,2))
    return np.array([np.diag([1,1,-1,-1]),*[np.block([[Z,s],[-s,Z]]) for s in sig]])

def spin_current_trace(f,fb):
    gam=make_gamma()
    slash=lambda p: np.einsum('i,ijk->jk',ETA@p,gam)
    return np.array([[np.trace(slash(fb)@g@slash(f)@h) for h in gam] for g in gam])

def numeric(seed,pairs):
    rng=np.random.default_rng(seed)
    xis=(0.37,1.,2.1,5.2)
    maxima={}; minimum={'propagator_gap': float('inf'),'root_offshell_gap':float('inf')}
    def observe(n,v): maxima[n]=max(maxima.get(n,0.),float(np.max(np.abs(v))))
    samples=[]; rejected=0
    for _ in range(pairs):
        while True:
            a=rng.uniform(.3,.8);m=rng.uniform(.7,1.35);M=rng.uniform(.4,1.1)
            b=rng.normal(size=3);b*=rng.uniform(.05,.55)/np.linalg.norm(b); B=boost(b)
            K=B@np.array([m,0.,0.,0.]); basis=(B@BASIS[1:].T).T
            n=rng.normal(size=3);n/=np.linalg.norm(n)
            f=B@(m/2*np.r_[1.,n]); fb=B@(m/2*np.r_[1.,-n])
            p,r=(np.r_[-np.sqrt(M*M+k@k),k] for k in rng.normal(size=(2,3)))
            q=-(K+p+r)
            Qs=[q+p,q+r]; H=p+r
            gaps=[abs(dot(Q,Q)-xi*m*m) for Q in Qs for xi in xis]+[abs(dot(H,H)-M*M)]
            off=abs(dot(q,q)-m*m)
            if min(gaps)>.08 and off>.2: break
            rejected+=1
        minimum['propagator_gap']=min(minimum['propagator_gap'],min(gaps))
        minimum['root_offshell_gap']=min(minimum['root_offshell_gap'],off)
        trace=spin_current_trace(f,fb)
        observe('mass_and_conservation', [dot(K,K)-m*m,dot(f,f),dot(fb,fb),dot(p,p)-M*M,dot(r,r)-M*M,*((K-f-fb))])
        observe('polarization_orthogonality', basis@ETA@K)
        observe('polarization_completeness', basis.T@basis-(-ETA+np.outer(K,K)/m**2))
        observe('left_current_conservation', trace@ETA@K)
        row=dict(g=2*a,m=m,M=M,xi=list(xis),q=q.tolist(),K=K.tolist(),p=p.tolist(),r=r.tolist(),fermion=f.tolist(),antifermion=fb.tolist(),checks=[])
        for xi in xis:
            diffs=[]; V=[]; G=[]
            for e in basis:
                numerical=direct_amplitudes(q,K,p,r,e,a,m,M,xi)
                ana=analytic_amplitudes(q,K,p,r,e,a,m,M,xi)
                scale=max(1.,np.max(np.abs(numerical)),np.max(np.abs(ana)))
                observe('per_graph_action_vs_closed_form_scaled',(numerical-ana)/scale)
                d=numerical[0]-1j*m*numerical[1]
                observe('contracted_st_scaled',d.sum()/scale)
                observe('channel_0_regular_scaled',(d[:2].sum()+2j*a*a*dot(q+p,e))/scale)
                observe('channel_1_regular_scaled',(d[2:4].sum()+2j*a*a*dot(q+r,e))/scale)
                observe('higgs_st_scaled',d[4]/scale)
                observe('omit_vector_quartic_control_scaled',d[:-1].sum()/scale)
                observe('omit_chi_exchange_from_goldstone_current_control_scaled',(d.sum()+1j*m*numerical[1,[1,3]].sum())/scale)
                V.append(numerical[0].sum());G.append(numerical[1].sum())
                diffs.append(d.sum())
            # f fbar spin sums contracted with an explicit sum over the 3 W polarizations.
            spin= basis@ETA@trace@ETA@basis.T
            abs2=lambda z: float(np.real(np.conjugate(z)@spin@z))
            vv=abs2(np.array(V)); gg=abs2(1j*m*np.array(G)); dd=abs2(np.array(diffs))
            observe('nonzero_spin_summed_vector_residue',vv)
            observe('spin_summed_identity_scaled',abs(vv-gg)/max(1.,vv,gg))
            observe('spin_summed_residual_scaled',abs(dd)/max(1.,vv,gg))
            row['checks'].append(dict(xi=xi, vector_contracted_pole_spin_norm2=vv,m_goldstone_pole_spin_norm2=gg,difference_pole_spin_norm2=dd))
        samples.append(row)
    controls={'pass_scaled':2e-10,'omission_min':1e-3,'positive_residue_min':1e-3}
    for k in ('omit_vector_quartic_control_scaled','omit_chi_exchange_from_goldstone_current_control_scaled'):
        assert maxima[k]>controls['omission_min'], maxima
    assert maxima['nonzero_spin_summed_vector_residue']>controls['positive_residue_min']
    assert max(val for key,val in maxima.items() if not key.startswith(('omit_','nonzero_')))<controls['pass_scaled'], maxima
    return dict(seed=seed,accepted=pairs,rejected=rejected,minimum=minimum,maxima=maxima,thresholds=controls,samples=samples)

if __name__=='__main__':
    cli=argparse.ArgumentParser()
    cli.add_argument('--seed',type=int,default=1209)
    cli.add_argument('--points',type=int,default=24)
    cli.add_argument('--output',type=Path,required=True)
    args=cli.parse_args()
    assert args.points>0
    report={'schema_version':1,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'symbolic':symbolic(),'numerical':numeric(args.seed,args.points),
            'versions':{'python':platform.python_version(),'numpy':np.__version__,'sympy':sp.__version__}}
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:report['numerical'][k] for k in ('accepted','minimum','maxima')},indent=2))

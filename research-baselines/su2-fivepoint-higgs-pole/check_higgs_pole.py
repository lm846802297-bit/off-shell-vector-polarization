#!/usr/bin/env python3
"""Physical hh scalar pole, analytic factorization vs complex contour extraction."""
import argparse
import hashlib
import itertools
import json
import platform
from pathlib import Path
import numpy as np
import sympy as sp

ETA=np.diag([1.,-1.,-1.,-1.])
E=np.eye(4)
S=(np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1,-1]))
Z=np.zeros((2,2))
GAMMA=np.array([np.diag([1,1,-1,-1]),*[np.block([[Z,i],[-i,Z]]) for i in S]])

def dot(a,b): return a@ETA@b

def slash(a): return np.einsum('i,ijk->jk',ETA@a,GAMMA)

def boost(b):
    g=(1-b@b)**-.5
    B=np.eye(4)
    B[0,0]=g; B[0,1:]=g*b; B[1:,0]=g*b
    B[1:,1:]+=g*g/(g+1)*np.outer(b,b)
    return B

def symbolic():
    a,m,B,xi,K2,qK,HK,qJ,HJ,KJ=sp.symbols('a m B xi K2 qK HK qJ HJ KJ')
    D=K2-m*m; C=K2-xi*m*m; rho=(1-xi)/C
    Rw=2*sp.I*a*m*B*(qJ-rho*qK*KJ)/D
    Rchi=a*B*(qJ-HJ-rho*(qK-HK)*KJ)/D
    difference=sp.factor((Rw-sp.I*m*Rchi).subs({HJ:-qJ-KJ,HK:-qK-K2}))
    assert sp.simplify(difference+sp.I*a*m*B*xi*KJ/C)==0
    assert sp.simplify(difference.subs(KJ,0))==0
    assert sp.simplify(Rw.subs(KJ,0)-2*sp.I*a*m*B*qJ/D)==0
    # Dirac projection, treated as noncommutative on-shell ideal.
    A,F=sp.symbols('A F',commutative=False)
    p=sp.expand(A*(A+F)*F)
    p=p.subs({A**2:0,F**2:0},simultaneous=True)
    assert p==0
    return dict(unprojected_delta='-i a m B xi slash(K)/(K^2-xi m^2)',
                projected_delta='0',gamma_projected_polynomial=str(p),
                symbolic_w=str(Rw),symbolic_chi=str(Rchi))

def bare_cubic(legs,a,m,M):
    """Numerical third multilinear coefficient of the original scalar action."""
    v=m/a
    if len(legs)!=3: raise ValueError(len(legs))
    assert np.max(np.abs(sum(x[1] for x in legs)))<2e-10
    f=0j
    for signs in itertools.product((-1,1),repeat=3):
        h=0j;c=0j; dh=np.zeros(4,complex); dc=dh.copy(); W=dh.copy()
        for s,(k,p,e) in zip(signs,legs):
            if k=='h': h+=s; dh-=1j*s*p
            elif k=='c': c+=s; dc-=1j*s*p
            elif k=='W': W+=s*e
            else: raise ValueError(k)
        # Two real components of D Phi. Square bilinearly, not Hermitian-conjugate
        # after analytic continuation of the external labels.
        u=dh-a*W*c
        w=dc+a*W*(v+h)
        t=h*h+2*v*h+c*c
        L=(dot(u,u)+dot(w,w))/2-M*M*t*t/(8*v*v)
        f+=np.prod(signs)*L
    return 1j*f/8

def deformation(z,M,BH,n,K):
    H=BH@np.array([np.sqrt(M*M+z+0j),0,0,0])
    im=1j*np.sqrt(M*M-(M*M+z)/4)
    ep=np.sqrt(M*M+z+0j)/2
    p=BH@np.r_[ep,im*n]
    r=BH@np.r_[ep,-im*n]
    q=-(K+H)
    return q,H,p,r

def closed(q,H,K,a,m,M,xi):
    ks=dot(K,K); v=m/a; B=3*M*M/v
    pol=ETA-(1-xi)*np.outer(K,K)/(ks-xi*m*m)
    W=2j*a*m*B*np.einsum('i,ijk->jk',ETA@q@pol@ETA,GAMMA)/(ks-m*m)
    C=a*B*np.einsum('i,ijk->jk',ETA@(q-H)@pol@ETA,GAMMA)/(ks-m*m)
    return np.array([W,C])

def finite_amplitude(z,args):
    K,BH,n,a,m,M,xi=args
    q,H,p,r=deformation(z,M,BH,n,K)
    assert max(abs(dot(p,p)-M*M),abs(dot(r,r)-M*M),abs(dot(H,H)-M*M-z))<2e-10
    ks=dot(K,K)
    prop=-1j*(ETA-(1-xi)*np.outer(K,K)/(ks-xi*m*m))/(ks-m*m)
    hhh=bare_cubic([('h',p,None),('h',r,None),('h',-H,None)],a,m,M)
    roots=[]
    for label,e in (('W',q),('c',None)):
        right=np.array([bare_cubic([(label,q,e),('W',K,u),('h',H,None)],a,m,M) for u in E])
        roots.append(np.einsum('i,ijk->jk',right@prop@ETA,1j*GAMMA)*(1j/z)*hhh)
    return np.array(roots)

def integrate(R,N,args):
    out=np.zeros((2,4,4),complex)
    # Use a phase-offset grid on each contour, neither contour includes the pole.
    for theta in 2*np.pi*(np.arange(N)+.371)/N:
        z=R*np.exp(1j*theta)
        out+=z*finite_amplitude(z,args)/N
    return out

def cm(a):
    return np.stack([np.real(a),np.imag(a)],axis=-1).tolist()

def evaluate(seed,count):
    rng=np.random.default_rng(seed)
    xis=[0.41,1.,2.35,4.8]
    maxmap={}; minimum={'K_propagator_gap':float('inf'),'root_offshell_gap':float('inf')}
    q2_range=[float('inf'),float('-inf')]
    def obs(name,x): maxmap[name]=max(maxmap.get(name,0.),float(np.max(np.abs(x))))
    archived=[]; rejects=0
    def randunit():
        x=rng.normal(size=3)
        return x/np.linalg.norm(x)
    for i in range(count):
        a=rng.uniform(.3,.8);m=rng.uniform(.7,1.45);M=rng.uniform(.55,1.5)
        while True:
            f=np.r_[rng.uniform(.45,2.),0.,0.,0.]
            f[1:]=f[0]*randunit()
            b=np.r_[rng.uniform(.45,2.),0.,0.,0.]
            b[1:]=b[0]*randunit()
            if i%2: b=-b  # crossed null leg tests spacelike K as well
            K=f+b
            BH=boost(rng.uniform(.03,.6)*randunit()); n=randunit()
            q,H,p,r=deformation(0,M,BH,n,K)
            gap=min(abs(dot(K,K)-x*m*m) for x in xis)
            off=abs(dot(q,q)-m*m)
            if gap>.18 and off>.15: break
            rejects+=1
        minimum['K_propagator_gap']=min(gap,minimum['K_propagator_gap'])
        minimum['root_offshell_gap']=min(off,minimum['root_offshell_gap'])
        q2=float(dot(q,q).real); q2_range=[min(q2_range[0],q2),max(q2_range[1],q2)]
        F=slash(f);A=slash(b)
        obs('external_mass_shell', [dot(f,f),dot(b,b),dot(p,p)-M*M,dot(r,r)-M*M,dot(H,H)-M*M,*((p+r-H))])
        obs('nonzero_scalar_complex_kinematics',max(abs(p.imag)))
        obs('projected_current_conserves',A@slash(K)@F)
        row=dict(a=a,m=m,M=M,q_squared=q2,crossed=bool(i%2),q=cm(q),H=cm(H),p=cm(p),r=cm(r),fermion=f.tolist(),antifermion=b.tolist(),gauges=[])
        ref=None
        for xi in xis:
            args=(K,BH,n,a,m,M,xi)
            first=integrate(.12*M*M,12,args)
            other=integrate(.055*M*M,20,args)
            expect=closed(q,H,K,a,m,M,xi)
            scale=max(1.,np.max(np.abs(expect)),np.max(np.abs(first)),np.max(np.abs(other)))
            d=first[0]-1j*m*first[1]
            dp=-1j*a*m*(3*M*M/(m/a))*xi*slash(K)/(dot(K,K)-xi*m*m)
            project=lambda x: np.array([A@y@F for y in x])
            P=project(first)
            Pscale=max(1.,np.max(np.abs(P)))
            if ref is None: ref=P
            for rho in [.12,.055]:
                z=rho*M*M*np.exp(.369j)
                _,Hd,pd,rd=deformation(z,M,BH,n,K)
                obs('deformed_mass_shell_and_pole', [dot(pd,pd)-M*M,dot(rd,rd)-M*M,dot(Hd,Hd)-M*M-z])
            obs('direct_contour_vs_factorization_scaled',(first-expect)/scale)
            obs('different_contour_radius_and_order_scaled',(first-other)/scale)
            obs('unprojected_delta_vs_slashK_scaled',(d-dp)/scale)
            obs('nonzero_unprojected_delta_scaled',d/scale)
            obs('projected_delta_scaled',A@d@F/Pscale)
            obs('projected_gauge_variation_scaled',(P-ref)/max(Pscale,np.max(np.abs(ref))))
            obs('nonzero_projected_signal',P)
            obs('wrong_phase_control_scaled',(P[0]+1j*m*P[1])/Pscale)
            rec=dict(xi=xi,nonprojected_st_max=float(np.max(np.abs(d))),
                     projected_st_max=float(np.max(np.abs(A@d@F))),projected_residue_max=float(np.max(np.abs(P))))
            if i==0:
                rec.update(raw_matrices_contour=cm(first),raw_matrices_factorized=cm(expect),projected_matrices=cm(P))
            row['gauges'].append(rec)
        archived.append(row)
    controls={'upper_scaled':2e-9,'strong_lower':1e-3}
    for name in ('nonzero_unprojected_delta_scaled','nonzero_projected_signal','nonzero_scalar_complex_kinematics','wrong_phase_control_scaled'):
        assert maxmap[name]>controls['strong_lower'],maxmap
    assert max(v for k,v in maxmap.items() if not k.startswith(('nonzero_','wrong_phase_')))<controls['upper_scaled'],maxmap
    return dict(seed=seed,accepted=count,rejected=rejects,q_squared_range=q2_range,minimum=minimum,maxima=maxmap,thresholds=controls,samples=archived)

def main():
    cli=argparse.ArgumentParser()
    cli.add_argument('--seed',type=int,default=20261123)
    cli.add_argument('--points',type=int,default=16)
    cli.add_argument('--output',type=Path,required=True)
    a=cli.parse_args()
    assert a.points>0
    data=dict(schema=1,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),symbolic=symbolic(),
              numerical=evaluate(a.seed,a.points),versions=dict(python=platform.python_version(),numpy=np.__version__,sympy=sp.__version__))
    a.output.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:data['numerical'][k] for k in ('accepted','rejected','q_squared_range','minimum','maxima')},indent=2))
if __name__=='__main__':
    main()

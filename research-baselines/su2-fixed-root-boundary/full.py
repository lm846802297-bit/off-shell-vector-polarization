import hashlib,sys,json
from pathlib import Path
import numpy as np
from boundary import closed
import importlib.util
LOCAL=Path(__file__).resolve().parent
PARENT=Path(sys.argv[1]).resolve()
PIN='196974fe69fa540d3b3b364cda2f5d11933ca79a12616188401a36c999a05e0e'
path=PARENT/'su2-sixpoint-connected-tree'/'check_connected_tree.py'
assert hashlib.sha256(path.read_bytes()).hexdigest()==PIN
spec=importlib.util.spec_from_file_location('parenttree',path)
orig=importlib.util.module_from_spec(spec)
spec.loader.exec_module(orig)
dot,slash,GAM,PL,I=orig.dot,orig.slash,orig.GAM,orig.PL,orig.I
cm=orig.fixture.cm
def spin(vb,G,u):return np.array([vb@gg@G@u for gg in GAM])

def decomp(params,data,e,n):
    a,m,H,mu,q,K,p,r,_=params
    f,F,u,vb=data
    b=mu*a/m;v=m/a
    j=spin(vb,PL,u)
    end=-a*e/(dot(K,K)-m*m)
    init=-a*(j-q*dot(q,j)/(m*m))/(dot(q,q)-m*m)
    wh=2j*a*m
    B=3*H*H/(v*(dot(p+r,p+r)-H*H))
    S=lambda l:(slash(l)+mu*I)/(dot(l,l)-mu*mu)
    def hard(P,h):
        L=2*dot(P,h);D=dot(P,P)-mu*mu
        return slash(h)/L,(slash(P)+mu*I)/L-slash(h)*D/L**2
    bridge=1j*closed(q,K,p,r,init,end,n,a,m,H)[0]
    mix=0j
    for x in (p,r):
        k=q+x
        s0,s1=hard(f+x,n)
        t0,t1=hard(f-k,-n)
        J0=np.array([vb@(g@PL@s0+t0@g@PL)@u for g in GAM])
        J1=np.array([vb@(g@PL@s1+t1@g@PL)@u for g in GAM])
        assert abs(dot(n,J0))<1e-11*max(1,np.linalg.norm(J0)*np.linalg.norm(n))
        mix+=wh*a*b*dot(end,n)*(dot(k,J0)+dot(n,J1))/(m*m*2*dot(k,n))
    U0=hard(f+p,n)[0]
    V0=hard(f+K+p,n)[0]
    W=slash(end)@PL
    double=-2j*a*b*b*(vb@(W@S(f+p+r)@U0+V0@W@U0+V0@S(f+K)@W)@u)
    split=-1j*a*b*B*(vb@(W@S(f+p+r)+S(f+K)@W)@u)
    return dict(bridge=bridge,mixed=mix,double=double,split=split)

def moving_poles(graphs,P,shift,m,H,mu,all_gauges):
    roots=[]
    for t in graphs:
        for _,kids in orig.walk(t):
            for k,s in kids:
                if isinstance(s,int):continue
                leaf=list(orig.leaves(s)); l=P[leaf].sum(axis=0); n=shift[leaf].sum(axis=0)
                if np.linalg.norm(n)<1e-12:continue
                D=2*dot(l,n)
                assert abs(D)>1e-5
                masses=((m*m,.63*m*m,1.7*m*m) if all_gauges else (m*m,)) if k in 'wc' else [H*H] if k=='h' else [mu*mu] if k in 'fF' else [0.]
                roots += [(mm-dot(l,l))/D for mm in masses]
    return list(dict.fromkeys(roots))

def cv(a):
    ar=np.asarray(a)
    return ar[...,0]+1j*ar[...,1]

def main():
    graphs=orig.gtree();kinds=np.array([orig.kind(t) for t in graphs])
    rows=[]
    for label in ('primary','holdout'):
        src=PARENT/'su2-sixpoint-connected-tree'/f'results-{label}.json'
        for ii,ex in enumerate(json.loads(src.read_text())['fixtures']):
            a,m,H,mu=[ex['parameters'][k] for k in ('a','m','M','mu')]
            P=cv(ex['external_momenta']).real
            f,F,g,G,p,r=P
            dat=cv(ex['external_spinors']);u,vb,ug,vbg=dat
            q=f+F;K=g+G;e=spin(vbg,I,ug)
            hh=p+r
            boost=orig.fixture.boost(hh[1:]/hh[0])
            pin=orig.fixture.boost(-hh[1:]/hh[0])@p
            k=pin[1:]/np.linalg.norm(pin[1:])
            ax=np.eye(3)[np.argmin(abs(k))]; orth=ax-(k@ax)*k;orth/=np.linalg.norm(orth)
            B0={k:cv(v) for k,v in ex['unitary_classes'].items()}
            at0=complex(sum(B0.values()))
            for sig in (1,-1):
                sh=boost@np.r_[0.,orth+1j*sig*np.cross(k,orth)]
                assert max(abs(dot(sh,p)),abs(dot(sh,r)),abs(dot(sh,sh)))<1e-12
                delta=np.zeros_like(P,dtype=complex);delta[4]=sh;delta[5]=-sh
                vals=decomp((a,m,H,mu,q,K,p,r,0), (f,F,u,vb),e,sh)
                roots=moving_poles(graphs,P,delta,m,H,mu,False)
                roots=[r0 for i,r0 in enumerate(roots) if all(abs(r0-r1)>1e-9 for r1 in roots[:i])]
                rad=3.8*max([1.]+list(map(abs,moving_poles(graphs,P,delta,m,H,mu,True))))
                def unitary(z):
                    return orig.unitary_closed((a,m,H,mu,q,K,p+z*sh,r-z*sh,0),(f,F,u,vb),e)
                pole_sum=0j
                for i,r0 in enumerate(roots):
                    small=.12*min([abs(r0)]+[abs(r0-r1) for j,r1 in enumerate(roots) if i!=j])
                    assert small>1e-7
                    ring=small*np.exp(2j*np.pi*np.arange(24)/24)
                    pole_sum-=np.mean([dz*sum(unitary(r0+dz).values()) for dz in ring])/r0

                node=rad*np.exp(2j*np.pi*np.arange(48)/48)
                ud={k:[] for k in vals}
                for z in node:
                    uz=unitary(z)
                    for key,x in uz.items():ud[key].append(x)
                bdiff={k:abs(vals[k]-np.mean(ud[k]))/max(1.,abs(vals[k])) for k in vals}
                assert max(bdiff.values())<1e-8,(label,ii,sig,bdiff)
                bsum=sum(vals.values())
                close=abs(at0-pole_sum-bsum)/max(1.,abs(at0),abs(bsum))
                assert close<1e-8,(label,ii,sig,close)

                extras={}
                if ii in (0,1):
                    ext=tuple([*dat,np.ones(1),np.ones(1)])
                    for xi in (.63,1.,1.7):
                        direct=np.array([orig.expanded_sum(graphs,P+z*delta,ext,a,m,H,mu,xi).sum() for z in node[::3]])
                        cvxi=abs(bsum-direct.mean())/max(1.,abs(bsum))
                        np0=orig.expanded_sum(graphs,P,ext,a,m,H,mu,xi).sum()
                        cv0=abs(np0-at0)/max(1.,abs(at0))
                        assert max(cvxi,cv0)<1e-8,(label,ii,sig,xi,cvxi,cv0)
                        extras[str(xi)]=[cvxi,cv0]
                rows.append(dict(input=label,index=ii,sign=sig,radius=rad,classes={k:cm(v) for k,v in vals.items()},
                     ground_amplitude=cm(at0),boundary=cm(bsum),physical_pole_sum=cm(pole_sum),moving_poles=len(roots),reconstruction_residual=close,class_boundary_residual=bdiff,Rxi_outer_and_zero=extras))
    print(json.dumps(dict(parent_commit='ba7dfb747da24809d1c1d9250bd9f8f6ce3677bf',cases=rows),indent=2))
if __name__=='__main__': main()

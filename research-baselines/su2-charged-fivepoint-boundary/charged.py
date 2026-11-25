import numpy as np
import json
from collections import Counter
from itertools import combinations,permutations
from functools import lru_cache
m=1.; M=.79; a=.31; v=m/a; gh=2*a*m
eta=np.diag([1.,-1.,-1.,-1.]); cyc=((0,1,2),(1,2,0),(2,0,1))
def d(x,y):return x@eta@y
def delta(x):return d(x,x)-m*m
def U(q,e):return (e-q*d(q,e)/m**2)/delta(q)
def V(p,e):
    return sum(d(e[i],e[j])*d(p[i]-p[j],e[k]) for i,j,k in cyc)
def c0(h):return 2*a*a*(1+3*M*M/(d(h,h)-M*M))
def direct(k,e,p,r):
    H=p+r; cs=c0(H); base=0j; join=0j
    for i,j,s in cyc:
        K=k[i]+H
        same=U(K,cs*e[i]+gh**2*(U(k[i]+p,e[i])+U(k[i]+r,e[i])))
        base+=V((K,k[j],k[s]),(same,e[j],e[s]))
        for x,y in ((p,r),(r,p)):
            l,jl=k[i]+x,k[j]+y
            join+=gh**2*V((l,jl,k[s]),(U(l,e[i]),U(jl,e[j]),e[s]))
    return np.array([base,join])

def boundary(k,e,p,r,n):
    H=p+r; cs=c0(H); t0=0j; bridge=0j; cross=0j
    ep=[d(n,u) for u in e]; kl=[2*d(n,q) for q in k]
    assert min(map(abs,kl))>1e-8
    R=np.array(ep)/np.array(kl)
    lead=0j
    for i,j,s in cyc:
        S=np.zeros(4,dtype=complex)
        for x in (p,r):
            l=k[i]+x
            S-= (l*ep[i]+n*d(l,e[i])-n*R[i]*delta(l))/(m**2*kl[i])
        K=k[i]+H
        t0+=V((K,k[j],k[s]),(U(K,cs*e[i]),e[j],e[s]))
        bridge+=gh**2*V((K,k[j],k[s]),(U(K,S),e[j],e[s]))
        for sign,x,y in ((1,p,r),(-1,r,p)):
            l,jl=k[i]+x,k[j]+y
            di=(d(l,e[i])-R[i]*delta(l))/kl[i]
            dj=-(d(jl,e[j])-R[j]*delta(jl))/kl[j]
            cross+=gh**2/m**2*(R[i]*d(e[j],e[s])-R[j]*d(e[i],e[s])
                   +R[i]*R[j]*d(l,e[s])+d(n,e[s])*(R[i]*dj+di*R[j]))
            lead+=sign*gh**2/m**2*R[i]*R[j]*d(n,e[s])
    assert abs(lead)<1e-11
    return np.array([t0,bridge,cross])

def cross_compact(k,e,p,r,n):
    H=p+r
    slopes=np.array([2*d(q,n) for q in k])
    en=np.array([d(x,n) for x in e])
    R=en/slopes
    F=np.array([(d(H,u)-Ri*(2*M*M+2*d(q,H)))/si for q,u,Ri,si in zip(k,e,R,slopes)])
    return gh**2/m**2*sum(R[i]*R[j]*d(k[i]-k[j],e[s])+
         en[s]*(R[j]*F[i]-R[i]*F[j]) for i,j,s in cyc)

def physical_basis(q):
    out=[]
    for u0 in np.eye(4,dtype=complex):
        u=u0-q*d(q,u0)/m**2
        for z in out:u=u+z*d(z,u)
        norm=d(u,u)
        if abs(norm)>1e-7:out.append(u/np.sqrt(-norm))
        if len(out)==3:break
    assert len(out)==3
    assert max(abs(d(q,u)) for u in out)<1e-9
    return out

def factored(k,e,x,y,i,n,sign):
    # The first subprocess is W_i h -> W. On the other side the
    # remaining Higgs is attached to any of the three arms of the cubic vertex.
    L=k[i]+x
    ks=np.array(k,dtype=complex);ks[i]=L
    ans=0j
    for eps in physical_basis(L):
        es=np.array(e,dtype=complex);es[i]=eps
        a4=0j
        for j in range(3):
            K=ks.copy();E=es.copy()
            K[j]=ks[j]+y;E[j]=gh*U(K[j],es[j])
            a4+=V(K,E)
        ans-=gh*d(e[i],eps)*a4
    return ans/(sign*2*d(k[i],n))

def n_orthogonal(k,e,n,count):
    E=e.copy()
    for i in range(count):
        mat=np.array([k[i]@eta,n.real@eta,n.imag@eta])
        _,_,vv=np.linalg.svd(mat)
        E[i]=vv[-1]/np.sqrt(-d(vv[-1],vv[-1]))
    return E

def boost(b):
    G=1/np.sqrt(1-b@b); A=np.eye(4)
    A[0,0]=G;A[0,1:]=A[1:,0]=G*b
    A[1:,1:]+=G*G/(G+1)*np.outer(b,b)
    return A
def sample(rng):
    u,w=rng.normal(size=(2,3)); t=-(u+w)
    eu,ew,et=[np.sqrt(h*h+q@q) for q,h in [(u,M),(w,M),(t,m)]]
    E=(eu+ew+et)/2; cx=np.sqrt(E*E-m*m)
    k=np.array([[E,0,0,cx],[E,0,0,-cx],[-et,*(-t)]])
    p=-np.array([eu,*u]);r=-np.array([ew,*w])
    e=[]
    for x in k:
        a0=rng.normal(size=3)
        y=np.r_[x[1:]@a0/x[0],a0];e.append(y/np.sqrt(-d(y,y)))
    H=p+r; b=H[1:]/H[0]
    hp=boost(-b)@p; z=hp[1:]/np.linalg.norm(hp[1:])
    x=np.eye(3)[np.argmin(abs(z))];x=x-z*(z@x);x/=np.linalg.norm(x)
    y=np.cross(z,x)
    ns=[boost(b)@np.r_[0,x+s*1j*y] for s in (1,-1)]
    return k,np.array(e),p,r,ns

# Unordered labelled partitions and independent rooted color/type tree census.
def part(idx,n):
    if n==1:
        yield (idx,)
        return
    for j in range(1,len(idx)-n+2):
        for group in combinations(idx[1:],j-1):
            first=(idx[0],)+group
            for tail in part(tuple(x for x in idx if x not in first),n-1):
                yield (first,)+tail

TYPES=('1','2','3','h','h')
# Field-content rules, including vertices that could have contributed.
RULES=[('1','2','3'),('h','h','h'),('h','h','h','h')]
RULES += [(w,w,'h') for w in '123']+[(w,w,'h','h') for w in '123']
RULES += [(u,u,w,w) for u,w in combinations('123',2)]
@lru_cache(None)
def trees(idx,w):
    if len(idx)==1:
        return (idx[0],) if TYPES[idx[0]]==w else ()
    result=[]
    import itertools
    for v in RULES:
        if w not in v:continue
        inputs=list(v);inputs.remove(w)
        for p in part(idx,len(inputs)):
            for ts in sorted(set(permutations(inputs))):
                for children in itertools.product(*(trees(i,t) for i,t in zip(p,ts))):
                    result.append((w,tuple(zip(ts,children))))
    return tuple(result)
def graph_class(T):
    types=[]
    def walk(T):
        if isinstance(T,int):return
        w,kids=T
        types.append(tuple(sorted((w,*(a for a,b in kids)))))
        for _,q in kids:walk(q)
    walk(T)
    assert types.count(('1','2','3'))==1
    return 'split' if ('h','h','h') in types else 'seagull' if any(len(t)==4 for t in types) else 'ordered'

def raw(graphs,k,e,p,r):
    # Independent rooted graph evaluation. Cubic gauge rule is the cyclic
    # covector tensor generated separately from the three derivative slots.
    P=[*k,p,r]; E=[*e,1.,1.]; basis=np.eye(4)
    def prop(w,l,y):
        if w=='h':return y/(d(l,l)-M*M)
        return (y-l*d(l,y)/m**2)/delta(l)
    def accum(t):
        if isinstance(t,int):return P[t],E[t]
        w,kids=t; q=[]; Vext=[]
        for wi,ti in kids:
            l,f=accum(ti)
            q.append(l)
            Vext.append(f if isinstance(ti,int) else prop(wi,l,f))
        x=sum(q)
        ts=(w,*(a for a,b in kids))
        if set(ts)=={'1','2','3'}:
            def term(a0,a1,a2):
                ai=[a0,a1,a2]; mom=[-x,*q]
                total=0
                for u,j,s in permutations(range(3)):
                    c=tuple(int(ts[t])-1 for t in (u,j,s))
                    f=1 if c in cyc else -1
                    total+=f*d(ai[u],ai[j])*d(mom[u],ai[s])
                return total
            return x,eta@np.array([term(y,*Vext) for y in basis])
        if all(t=='h' for t in ts):
            return x,3*M*M/v*np.prod(Vext)
        if w=='h':
            hv=np.prod([y for (_t,_i),y in zip(kids,Vext) if _t=='h'])
            vert=[y for (_t,_i),y in zip(kids,Vext) if _t!='h']
            return x,(gh if len(ts)==3 else 2*a*a)*d(*vert)*hv
        hv=np.prod([y for (_t,_i),y in zip(kids,Vext) if _t=='h'])
        vert=[y for (_t,_i),y in zip(kids,Vext) if _t!='h']
        return x,(gh if len(ts)==3 else 2*a*a)*vert[0]*hv
    return sum(d(e[2],accum(t)[1]) for t in graphs)

def circle(fun,center,R,N=64,small=False):
    z=R*np.exp(2j*np.pi*np.arange(N)/N)
    return np.mean([fun(center+t)*(t if small else 1) for t in z],axis=0)
def allroots(k,p,r,n):
    return [-delta(q+x)/(sign*2*d(q,n)) for q in k for x,sign in ((p,1),(r,-1))]
def check_control(k,e,p,r,n):
    # B is linear in e3. Choose physical complex polarization in its kernel.
    M0=np.array([k[2]@eta])
    _,_,V=np.linalg.svd(M0)
    basis=V[1:]
    coeff=np.array([boundary(k,np.array([e[0],e[1],z]),p,r,n).sum() for z in basis])
    A=np.array([direct(k,np.array([e[0],e[1],z]),p,r).sum() for z in basis])
    _,_,N=np.linalg.svd(coeff[None,:])
    roots=N[1:].conj()
    c=roots[np.argmax(abs(roots@A))]
    ee=np.array([e[0],e[1],c@basis]); B=boundary(k,ee,p,r,n).sum(); az=direct(k,ee,p,r).sum()
    assert abs(B)<1e-10 and abs(az)>1e-6,(B,az)
    return ee,B,az

def cv(a):
    return np.stack((np.real(a),np.imag(a)),axis=-1).tolist()
def main():
    graphs=trees((0,1,3,4),'3')
    cnt=Counter(map(graph_class,graphs))
    assert len(set(graphs))==18 and cnt==dict(ordered=12,seagull=3,split=3)
    outs=[]
    for seed in (20261125,171108):
        rng=np.random.default_rng(seed)
        for j in range(10):
            k,e,p,r,sh=sample(rng)
            assert np.linalg.norm(k.sum(axis=0)+p+r)<1e-10
            assert max(abs(d(x,y)) for x,y in zip(k,e))<1e-10
            for sign,n in zip((1,-1),sh):
                assert max(abs(d(n,p)),abs(d(n,r)),abs(d(n,n)))<1e-10
                Rts=allroots(k,p,r,n)
                diffs=[abs(x-y) for i,x in enumerate(Rts) for y in Rts[:i]]
                assert min(diffs)>1e-3
                R=4*max([1.]+list(map(abs,Rts)))
                def fn(z,ee=e):return direct(k,ee,p+z*n,r-z*n)
                B=boundary(k,e,p,r,n); A=fn(0.).sum()
                cC=circle(fn,0,R)
                Cbar=np.array([B[0]+B[1],B[2]])
                sc=max(1.,abs(A),*abs(B))
                br=max(abs(cC-Cbar))/sc
                assert br<1e-9
                SS=0j; factor=0.
                for i,z in enumerate(Rts):
                    sm=.08*min([abs(z)]+[abs(z-y) for l,y in enumerate(Rts) if l!=i])
                    rho=circle(fn,z,sm,40,True).sum()
                    SS-=rho/z
                    sign_h=1 if i%2==0 else -1
                    ph,rr=p+z*n,r-z*n
                    if sign_h==-1:ph,rr=rr,ph
                    fact=factored(k,e,ph,rr,i//2,n,sign_h)
                    factor=max(factor,abs(fact-rho)/max(1.,abs(fact),abs(rho)))
                assert factor<1e-9,(seed,j,factor)
                compact=abs(cross_compact(k,e,p,r,n)-B[2])/sc
                assert compact<1e-12
                cr=abs(A-SS-B.sum())/sc
                assert cr<1e-9,(seed,j,cr)
                # No growing coefficients O(z), O(z^2), O(z^3) in the complete amplitude.
                grow=max(abs(circle(lambda z:fn(z).sum()/(z/R)**q,0,R)) for q in (1,2,3))/sc
                assert grow<1e-9
                gr=max(abs(raw(graphs,k,e,p+z*n,r-z*n)-fn(z).sum()) for z in (0.,.19+.72j))/(sc)
                assert gr<1e-9,(seed,j,gr,raw(graphs,k,e,p,r),A)
                row=dict(seed=seed,sample=j,sign=sign,k=cv(k),e=cv(e),h=cv([p,r]),n=cv(n),terms=cv(B),A0=cv(A),poles_only=cv(SS),
                         radius=R,roots=cv(Rts),scaled_cn_err=br,scaled_recon_err=cr,scaled_grow=grow,scaled_graph_oracle=gr,scaled_factorization=factor,scaled_compact=compact)
                if j==0:
                    e_one=n_orthogonal(k,e,n,1);e_two=n_orthogonal(k,e,n,2)
                    b_one=boundary(k,e_one,p,r,n);b_two=boundary(k,e_two,p,r,n)
                    assert abs(b_one[2])>1e-6 and abs(b_two[2])<1e-10 and abs(b_two.sum())>1e-6
                    row['one_two_n_orthogonal_control']=dict(one_cross=cv(b_one[2]),two_cross=cv(b_two[2]),two_total=cv(b_two.sum()))
                    ee,bb,az=check_control(k,e,p,r,n)
                    assert abs(circle(lambda z:fn(z,ee).sum(),0,R))<1e-9
                    row['nontrivial_zero_boundary_control']=dict(e=cv(ee),boundary=cv(bb),A0=cv(az))
                outs.append(row)
    print(json.dumps(dict(diagrams=dict(cnt),m=m,M=M,a=a,rows=outs),separators=(',',':')))
if __name__=='__main__':main()

"""Independent labelled finite-Rxi trees from the Cartesian SU(2) action."""
import numpy as np
from itertools import combinations, permutations, product
from collections import Counter
from functools import lru_cache
import json
import charged as b
m,M,a,v=b.m,b.M,b.a,b.v
g,gh=2*a,b.gh
eta=np.array([1.,-1.,-1.,-1.])
eye=np.eye(4)
cyc=((0,1,2),(1,2,0),(2,0,1))
W=tuple('ABC'); G=tuple('xyz')
color={t:j for fs in (W,G) for j,t in enumerate(fs)}
TYPES=('A','B','C','h','h')
def d(u,w): return np.sum(u*eta*w)
def eps(i,j,k): return 0 if len({i,j,k})!=3 else 1 if (i,j,k) in cyc else -1
# A catalogue, not the unitary-gauge topology list. Include inactive vertices.
V=[('WWW',W),('hhh',('h',)*3),('hhhh',('h',)*4)]
for i in range(3):
    w,q=W[i],G[i]
    V.extend([('WWh',(w,w,'h')),('WWhh',(w,w,'h','h')),('WhG',(w,'h',q)),
        ('WGG',(w,G[(i+1)%3],G[(i+2)%3])),('hGG',('h',q,q)),
        ('hhGG',('h','h',q,q)),('GGGG',(q,)*4)])
    V += [('WWGG',(w,w,t,t)) for t in G]
for i,j in combinations(range(3),2):
    V.extend([('WWWW',(W[i],W[i],W[j],W[j])),('GGGG',(G[i],G[i],G[j],G[j]))])

def split(items,n):
    if n==1:
        yield (items,); return
    # root by smallest label in the first block to avoid overcounting
    for size in range(1,len(items)-n+2):
        for X in combinations(items[1:], size-1):
            A=(items[0],*X)
            for B in split(tuple(i for i in items if i not in A),n-1):
                yield (A,*B)
@lru_cache(None)
def build(labels,target):
    if len(labels)==1: return (labels[0],) if TYPES[labels[0]]==target else ()
    all=[]
    for tag,types in V:
        if target not in types:continue
        leaves=list(types);leaves.remove(target)
        for p in split(labels,len(leaves)):
            for fs in sorted(set(permutations(leaves))):
                for ch in product(*(build(X,t) for X,t in zip(p,fs))):
                    all.append((target,tag,tuple(zip(fs,ch))))
    return tuple(all)
GRAPHS=build((0,1,3,4),'C')

def vertex(tag,f,p,e):
    if tag=='WWW':
        return g*sum(eps(color[f[i]],color[f[j]],color[f[k]])*d(e[i],e[j])*d(p[i]-p[j],e[k]) for i,j,k in cyc)
    wi=[i for i,t in enumerate(f) if t in W]
    scal=np.prod([e[j] for j in range(len(f)) if j not in wi])
    if tag in ('WWh','WWhh','WWGG'): return 1j*(gh if tag=='WWh' else 2*a*a)*d(*(e[i] for i in wi))*scal
    if tag=='WhG':
        j=f.index('h');i=wi[0];k=next(l for l,t in enumerate(f) if t in G)
        return a*d(p[j]-p[k],e[i])*scal
    if tag=='WGG':
        i=wi[0];j,k=[l for l in range(3) if l!=i]
        return a*eps(color[f[i]],color[f[j]],color[f[k]])*d(p[k]-p[j],e[i])*scal
    if tag in ('hhh','hGG','hhhh','hhGG','GGGG'):
        n=3 if tag in ('hhh','hhhh') or (tag=='GGGG' and len(set(f))==1) else 1
        return -1j*n*M*M/(v if len(f)==3 else v*v)*scal
    if tag=='WWWW':
        return -1j*g*g/4*sum(sum(eps(color[f[i]],color[f[j]],z)*eps(color[f[k]],color[f[l]],z) for z in range(3))*d(e[i],e[k])*d(e[j],e[l]) for i,j,k,l in permutations(range(4)))
    raise ValueError(tag)
def prop(f,p,x,xi):
    p2=d(p,p)
    if f in W:
        # use the standard Rxi form directly, not a unitary subtraction
        return -1j*(x-(1-xi)*p*d(p,x)/(p2-xi*m*m))/(p2-m*m)
    return 1j*x/(p2-(M*M if f=='h' else xi*m*m))

def tags(t):
    if isinstance(t,int): return []
    return [t[1]]+[y for _,c in t[2] for y in tags(c)]
CLASSES=tuple(('doubleG' if 'WGG' in tags(t) else 'singleG' if 'WhG' in tags(t) else 'vector') for t in GRAPHS)

def eval_amp(k,e,p,r,xi):
    ps=[*k,p,r]; es=[*e,1,1]
    def go(t,isroot=False):
        if isinstance(t,int):return ps[t],es[t]
        typ,tag,ch=t
        kids=[go(T) for _,T in ch]
        q=sum(p for p,_ in kids)
        F=(typ,*(f for f,_ in ch)); P=(-q,*(p for p,_ in kids)); E=tuple(x for _,x in kids)
        val=lambda u:vertex(tag,F,P,(u,*E))
        if isroot:return val(es[2])/g
        x=eta*np.array([val(u) for u in eye]) if typ in W else val(1)
        return q,prop(typ,q,x,xi)
    out=np.zeros(3,complex)
    for c,t in zip(CLASSES,GRAPHS): out[('vector','singleG','doubleG').index(c)]+=go(t,True)
    return out

def moment(fun,center,R,N,power=0):
    t=R*np.exp(2j*np.pi*np.arange(N)/N)
    return np.mean([fun(center+w)*w**power for w in t],axis=0)
def roots(k,p,r,n,m2):
    return [-(d(q+x,q+x)-m2)/(s*2*d(q,n))
            for q in k for x,s in ((p,1),(r,-1))]
def cv(x): return np.stack((np.real(x),np.imag(x)),axis=-1).tolist()
def main():
    assert len(set(GRAPHS))==30
    assert Counter(CLASSES)=={'vector':18,'singleG':6,'doubleG':6}
    rng=np.random.default_rng(2026112518); stats=[]; worst=dict(amplitude=0.,boundary=0.,spurious=0.,physical_residue=0.)
    seen=0.
    for index in range(6):
        k,e,p,r,nn=b.sample(rng)
        if index%2:
            q=k[0]; en=np.linalg.norm(q[1:]);e[0]=np.r_[en,q[0]*q[1:]/en]/m
        for side,n in enumerate(nn):
            pu=roots(k,p,r,n,m*m)
            c0=b.boundary(k,e,p,r,n).sum()
            for xi in (0.31,1.7,4.2):
                un=roots(k,p,r,n,xi*m*m); allp=pu+un
                R=2.7*max(1.,*map(abs,allp))
                fun=lambda z:eval_amp(k,e,p+z*n,r-z*n,xi)
                amps=[]
                for z in (0j,.13+.29j,-.74+.52j,1.4-.9j):
                    parts=fun(z); A=b.direct(k,e,p+z*n,r-z*n).sum()
                    err=abs(parts.sum()-A)/max(1.,abs(A),sum(abs(parts)))
                    worst['amplitude']=max(worst['amplitude'],err)
                    amps.append(dict(z=cv(z),parts=cv(parts),unitary=cv(A),absdiff=abs(parts.sum()-A)))
                bc=moment(fun,0,R,36)
                worst['boundary']=max(worst['boundary'],abs(bc.sum()-c0)/max(1.,abs(c0),sum(abs(bc))))
                res=[]
                for j,c in enumerate(allp):
                    radius=.17*min(1.,*(abs(c-d0) for i,d0 in enumerate(allp) if i!=j))
                    rr=moment(fun,c,radius,12,1)
                    if j<6:
                        ii,ipt=divmod(j,2); sg=1 if ipt==0 else -1
                        x,y=(p,r) if sg==1 else (r,p)
                        expected=b.factored(k,e,x+c*sg*n,y-c*sg*n,ii,n,sg)
                        key='physical_residue'
                    else:
                        expected=0.;key='spurious'
                        seen=max(seen,max(abs(rr)))
                    worst[key]=max(worst[key],abs(rr.sum()-expected)/max(1.,abs(expected),sum(abs(rr))))
                    res.append(dict(pole=cv(c),from_groups=cv(rr),expected=cv(expected),category=key))
                stats.append(dict(index=index,branch=side,xi=xi,k=cv(k),e=cv(e),p=cv(p),r=cv(r),n=cv(n),
                    amplitudes=amps,contour_radius=R,boundary=cv(bc),unitary_boundary=cv(c0),residues=res))
    assert max(worst.values())<1e-8 and seen>1e-3,(worst,seen)
    print(json.dumps(dict(gauges=[.31,1.7,4.2],graphs=dict(Counter(CLASSES)),worst_scaled=worst,
              max_noncancelled_component_of_spurious_residue=seen,fixtures=stats),indent=2,allow_nan=False))
if __name__=='__main__':main()

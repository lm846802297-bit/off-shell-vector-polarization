"""Exact symbolic checks: CPython 3 + SymPy; project momenta all incoming."""
import json
import sympy as s
v=s.Matrix; I=s.I; z=s.Symbol('z')
metric=s.diag(1,-1,-1,-1)
def dot(a,b): return (a.T*metric*b)[0]
def sq(a): return s.expand(dot(a,a))
def D(a): return sq(a)-1
def U(q,e): return (e-q*dot(q,e))/D(q)
def bc(f): return s.cancel(f).subs(z,1/z).series(z,0,1).removeO().expand().coeff(z,0).simplify()
def zero(a): assert s.simplify(a)==0,a
def base(H,M): return s.Rational(1,2)*(1+3*M**2/(sq(H)-M**2))
# m=1, g2=1, gh=1; independent scalar numerator long division of each moving channel.
def from_four(k,u,w,p,r,M,n):
    H=p+r; b=2*dot(k,n)
    C0=base(H,M)*dot(u,w); cs=[]
    for x,sign in ((p,1),(r,-1)):
        L=k+x; X=dot(u,L);Y=dot(w,L);N=dot(u,n);P=dot(w,n)
        f=(dot(u,w)-(X+sign*z*N)*(Y+sign*z*P))/(D(L)+sign*z*b)
        c=(N*P*D(L)/b**2-(N*Y+P*X)/b)
        lead=-sign*N*P/b
        numerator=s.Poly(s.cancel((f-c-lead*z)*(D(L)+sign*z*b)),z)
        assert numerator.degree()<=0
        zero(bc(f)-c)
        cs.append((f,c))
    return C0,cs

def V(k,e):
    return sum(dot(e[i],e[j])*dot(k[i]-k[j],e[l]) for i,j,l in ((0,1,2),(1,2,0),(2,0,1)))
cyc=((0,1,2),(1,2,0),(2,0,1))
def raw_cross(k,e,p,r,n):
    y=0
    for i,j,l in cyc:
        for x,t in ((p+z*n,r-z*n),(r-z*n,p+z*n)):
            q=list(k); f=list(e)
            q[i]+=x; q[j]+=t
            f[i]=U(q[i],e[i]); f[j]=U(q[j],e[j])
            y+=V(q,f)
    return y
def shift_check(k,e,p,r,M,n):
    assert sum(k, s.zeros(4,1))+p+r==s.zeros(4,1)
    for ki,ei in zip(k,e):
        zero(sq(ki)-1); zero(dot(ki,ei))
    for x in (p,r): zero(sq(x)-M**2); zero(dot(x,n))
    zero(sq(n))
def main():
    Q=s.Rational
    p=v([-Q(5,3),Q(4,3),0,0]);r=v([-Q(5,3),-Q(4,3),0,0])
    n=v([0,0,1,I]); k=v([Q(5,3),0,0,Q(4,3)]);w=v([0,0,1,0])
    shift_check([k,v([Q(5,3),0,0,-Q(4,3)])],[w,w],p,r,1,n)
    fixed,terms=from_four(k,w,w,p,r,1,n)
    c_on=bc(fixed+sum(f for f,c in terms))
    zero(c_on-fixed-sum(c for f,c in terms))
    H=p+r;P=k+H; b=2*dot(k,n);R=dot(w,n)/b
    # Sew a true cubic current, projected conservation is tested independently.
    K=(v([3,2,2,0]),v([3,-2,2,0]),v([-s.sqrt(17),0,-4,0]))
    H=-sum(K,s.zeros(4,1)); p=H/2+v([0,Q(1,2),0,0]);r=H-p
    M=s.sqrt(sq(p)); E=v([0,0,0,1]); es=v([0,1,0,0])
    shift_check(K,[E,E,es],p,r,M,n)
    children=[]
    for i,j,l in cyc:
        P=K[i]+H; momenta=(P,K[j],K[l]); vecs=[None,(E,E,es)[j],(E,E,es)[l]]
        J=metric*v([V(momenta,[a0,*vecs[1:]]) for a0 in map(v,s.eye(4).tolist())])
        zero(dot(J,P))
        uu=(E,E,es)[i]
        fixed,t=from_four(K[i],uu,J,p,r,M,n)
        source_constants=[]
        for x,sign in ((p,1),(r,-1)):
            L=K[i]+x+sign*z*n
            source_constants.append(bc(dot(J,U(P,U(L,uu)))))
        zero((fixed+sum(c for f,c in t))/D(P) - fixed/D(P)-sum(source_constants))
        children.append(s.sstr(s.simplify((fixed+sum(c for f,c in t))/D(P))))
    out=[]
    for name,ee,want in [
        ('one_n_orthogonal_nonzero',[E,E,es],Q(1,4)),
        ('one_n_orthogonal_kernel',[E,E,v([4,0,s.sqrt(17),I*s.sqrt(17)])],0),
        ('two_n_orthogonal',[v([2,3,0,0]),E,es],0)
    ]:
        shift_check(K,ee,p,r,M,n)
        for nn in (n,s.conjugate(n)):
            slopes=[2*dot(q,nn) for q in K]; N=[dot(nn,x) for x in ee]
            R=[a/b for a,b in zip(N,slopes)]
            T=slopes[1]*K[0]-slopes[0]*K[1]
            reduced=-2*R[0]*R[1]*dot(ee[2],T)/slopes[2]
            zero(reduced-want)
            full=bc(raw_cross(K,ee,p,r,nn))
            zero(full-reduced)
            out.append([name,s.sstr(nn[3]),s.sstr(full)])
    print(json.dumps({'four_point_fixed':s.sstr(base(v([-Q(10,3),0,0,0]),1)*dot(w,w)),
    'four_point_ordered_constants':[s.sstr(c) for _,c in terms],'four_point_total':s.sstr(c_on),
    'same_arm_exact_checks':len(children),'cross_exact_controls':out},indent=2))
if __name__=='__main__': main()

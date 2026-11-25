from fractions import Fraction
import json
import numpy as np
import sympy as s
eta=np.diag([1.,-1.,-1.,-1.])
def dot(a,b):return a@eta@b
def at(q,K,p,r,j,e,a,m,H):
    def fun(z,n,xi):
        total=0j
        for Q in (q+p+z*n,q+r-z*n):
            D=dot(Q,Q)-m*m
            C=dot(Q,Q)-xi*m*m
            # Direct W propagator contraction and the two derivative W-chi-h vertices.
            ww=4*a*a*m*m*(dot(j,e)-(1-xi)*dot(j,Q)*dot(e,Q)/C)/D
            gg=-a*a*dot(j,2*Q-q)*dot(e,2*Q+K)/C
            total+=ww+gg
        return total +static(q,K,a,H)*dot(j,e)
    return fun

def static(q,K,a,H):
    return 2*a*a*(1+3*H*H/(dot(q+K,q+K)-H*H))

def closed(q,K,p,r,j,e,n,a,m,H):
    d=2*dot(q,n)
    Q=[q+p,q+r]; D=[dot(x,x)-m*m for x in Q]; Qs=q-K
    jn,en=dot(j,n),dot(e,n)
    B=static(q,K,a,H)*dot(j,e)-4*a*a*((dot(j,Qs)*en+jn*dot(e,Qs))/d-jn*en*sum(D)/d**2)
    R=[]
    for x,Dx,ds in zip(Q,D,[d,-d]):
        z=-Dx/ds; X=x+z*n*np.sign(1 if ds==d else -1)
        R.append((z,4*a*a*(m*m*dot(j,e)-dot(j,X)*dot(e,X))/ds))
    return B,R

def symtest():
    z,d,D1,D2,A1,A2,L1,L2,T=s.symbols('z d D1 D2 A1 A2 L1 L2 T',nonzero=True)
    f=(A1+L1*z+T*z**2)/(D1+d*z)+(A2-L2*z+T*z**2)/(D2-d*z)
    B=(L1+L2)/d-T*(D1+D2)/d**2
    assert s.simplify(s.limit(f,z,s.oo)-B)==0
    for dd,ll,AA,ds in [(D1,L1,A1,d),(D2,-L2,A2,-d)]:
        zp=-dd/ds
        f-= (AA+ll*zp+T*zp*zp)/(ds*(z-zp))
    assert s.cancel(f-B)==0
    q2,k2,qk,h=s.symbols('q2 k2 qk h',nonzero=True)
    assert s.factor(s.Matrix([[q2,qk,h],[qk,k2,-h],[h,-h,0]]).det())==-h*h*(q2+k2+2*qk)
    # Goldstone and non-Feynman propagator term give the same metric-independent polynomial.
    Q2,ms,xi=s.symbols('Q2 ms xi')
    assert s.cancel((ms*(1-xi)/(Q2-ms)+1)/(Q2-xi*ms)-1/(Q2-ms))==0

def ctrace(fun,R,c=0,n=96):
    circle=R*np.exp(2j*np.pi*np.arange(n)/n)
    return np.mean([fun(c+d)*d if c!=0 else fun(d) for d in circle])

def cases():
    symtest()
    rng=np.random.default_rng(261125)
    a,m,H=2/7,2/3,1.
    rows=[]
    ratios=[]
    for E in (1.25,1.4,1.6):
        p=np.array([E,0,0,np.sqrt(E*E-H*H)])
        r=p*np.array([1,1,1,-1])
        for q in (np.array([1.5,2/3, .2,2/7]),np.array([2.4,.27,.39,-.31])):
            K=-q-p-r
            for n in (np.array([0,1,1j,0]),np.array([0,1,-1j,0])):
                # Epsilon-vector dual to (q,K,n) via 3x3 minors. No complex conjugates.
                M=np.array([q,K,n])@eta
                t=np.array([(-1)**i*np.linalg.det(np.delete(M,i,axis=1)) for i in range(4)])
                assert max(abs(dot(t,q)),abs(dot(t,K)),abs(dot(t,n)))<1e-13
                assert abs(dot(t,t))>.1
                j,e=rng.normal(size=(2,4))
                j-=q*dot(j,q)/dot(q,q)
                e-=K*dot(e,K)/dot(K,K)
                for kind,L,R in [('generic',j,e),('witness',t,t)]:
                    amp=at(q,K,p,r,L,R,a,m,H)
                    bound,poles=closed(q,K,p,r,L,R,n,a,m,H)
                    assert min(abs(poles[0][0]),abs(poles[1][0]),abs(poles[0][0]-poles[1][0]))>.04
                    if kind=='witness': assert abs(bound-static(q,K,a,H)*dot(t,t))<1e-11 and abs(bound)>.01
                    for xi in (.6,1.,2.1):
                        orig=amp(0.,n,xi); eval=lambda z:amp(z,n,xi)
                        scale=max(1.,abs(orig),abs(bound))
                        resid_sum=sum(-rho/zp for zp,rho in poles)
                        assert abs(orig-bound-resid_sum)/scale<1e-11
                        radius=4*max(abs(poles[0][0]),abs(poles[1][0]),1.)
                        Rs=.08*min(abs(poles[0][0]),abs(poles[1][0]),abs(poles[0][0]-poles[1][0]))
                        Bnum=ctrace(eval,radius)
                        Bnum2=ctrace(eval,1.23*radius)
                        rr=max(abs(ctrace(eval,Rs,zp)-rho)/scale for zp,rho in poles)
                        # Check an unphysical Goldstone pole in each ordering separately
                        # at the full gauge-completed sewn level; neither gauge term alone is zero.
                        Csh=[(xi-1)*m*m/(2*dot(q,n)),-(xi-1)*m*m/(2*dot(q,n))]
                        gr=[]
                        if xi!=1:
                            for (zp,_),k in zip(poles,Csh):
                                zg=zp+k
                                dr=min(abs(zg-z0) for z0,_ in poles)
                                if dr>.002:
                                    gr.append(abs(ctrace(eval,.07*dr,zg))/scale)
                        worst=max(abs(Bnum-bound)/scale,abs(Bnum2-bound)/scale,rr,max(gr,default=0.))
                        assert worst<2e-10,(kind,xi,worst)
                        if kind=='witness': assert abs(orig-resid_sum)>.01
                        rows.append(dict(E=E,q=q.tolist(), null_sign=n[2].imag, kind=kind, xi=xi,
                            boundary=[bound.real,bound.imag], full=[orig.real,orig.imag],
                            physical_residue_residual=rr,gauge_pole_residual=max(gr,default=0.),
                            boundary_residual=abs(Bnum-bound)/scale))
    print(json.dumps(dict(cases=rows),indent=2))
if __name__=='__main__':cases()

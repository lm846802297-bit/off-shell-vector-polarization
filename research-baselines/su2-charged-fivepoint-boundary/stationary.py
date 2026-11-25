from charged import *
def case(s):
    # Four-point physical momenta with an h* of mass sqrt(s); split it
    # into two on-shell Higgs momenta with analytic complex relative momentum.
    E=2.
    b=np.sqrt(E*E-m*m)
    EW=(4*E*E+m*m-s)/(4*E)
    EH=2*E-EW
    w=np.sqrt(EW*EW-m*m)
    k=np.array([[E,0,0,b],[E,0,0,-b],[-EW,-w,0,0]])
    H=-sum(k)
    e=np.array([[0,0,1.,0],[0,0,1.,0],[0,0,0,1.]])
    L=boost(H[1:]/H[0])
    # Use a fixed generic axis, then the complex null direction
    # in its orthogonal transverse plane.
    u=np.array([1.,2.,2.])/3
    x=np.array([2.,-1.,0.])/np.sqrt(5.)
    y=np.cross(u,x)
    p=L@np.r_[-np.sqrt(s)/2,complex(s/4-M*M)**.5*u]
    r=H-p
    n=L@np.r_[0.,x+1j*y]
    assert max(abs(d(p,p)-M*M),abs(d(r,r)-M*M),abs(d(n,p)),abs(d(n,r)),abs(d(n,n)))<1e-12
    a4=gh*2*b*sum(1/delta(z+H) for z in k)
    return k,e,p,r,n,a4
def main():
    *_,a4=case(M*M)
    expected=(3*M*M/v)*a4
    assert abs(expected)>.005
    rows=[]
    for off in (1e-3,1e-4,1e-5,1e-6):
        k,e,p,r,n,_=case(M*M+off)
        B=boundary(k,e,p,r,n)
        rows.append(dict(offset=off,expected=expected,boundary_scaled=cv(off*B.sum()),
                         residual=abs(off*B.sum()-expected)))
    assert rows[-1]['residual']<rows[0]['residual']/100 and rows[-1]['residual']<1e-5
    print(json.dumps(rows,indent=2))
if __name__=='__main__':main()

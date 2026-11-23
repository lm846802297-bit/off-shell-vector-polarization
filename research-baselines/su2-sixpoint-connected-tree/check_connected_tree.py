#!/usr/bin/env python3
"""Labelled six-point tree census, action-tensor sewing and independent unitary sum."""
import argparse
from collections import Counter
from functools import cache
import hashlib
import importlib.util
from itertools import permutations, product
import json
from pathlib import Path

import numpy as np
import sympy as sp

BASE = Path(__file__).resolve().parents[1]
HELPER = BASE/'su2-sixpoint-double-higgs-radiation'/'check_double_higgs.py'
PIN = '832d7a89d99c77aef808948d68ac98aa989d820b62c588b606e61b1599f37f24'
if hashlib.sha256(HELPER.read_bytes()).hexdigest() != PIN:
    raise RuntimeError("pinned external spinor fixture source has changed")
spec = importlib.util.spec_from_file_location('sixpoint_fixtures', HELPER)
fixture = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fixture)
ETA, GAM, G5, PL, I, slash, dot = fixture.ETA, fixture.GAM, fixture.G5, fixture.PL, fixture.I4, fixture.slash, fixture.dot

LABELS = ('f', 'F', 'g', 'G', 'h', 'h')
DUAL = dict(f='F', F='f', g='G', G='g', w='w', c='c', h='h')
VERTICES = tuple(tuple(sorted(v)) for v in (
    'Ffw', 'Ffc', 'Ffh', 'Ggw', 'wwh', 'wch', 'cch', 'hhh',
    'wwhh', 'wwcc', 'hhhh', 'hhcc', 'cccc'))
EXPECT = dict(bridge=11, mixed=8, double=6, split=2)
GAUGES = (.43, .71, 1., 1.8, 2.65)
ROOT = 3
EOS = 1e-9
NEG = 1e-4
GAP = .22


@cache
def parts(mask, n):
    """Unordered disjoint partitions, with canonical first set containing low bit."""
    if n == 0:
        return ((),) if mask == 0 else ()
    if n == 1:
        return ((mask,),) if mask else ()
    ans = []
    low = mask & -mask
    rem = mask ^ low
    q = rem
    while q:
        q = (q-1) & rem
        first = q | low
        ans.extend((first, *p) for p in parts(mask ^ first, n-1))
    return tuple(ans)


@cache
def trees(mask, incoming):
    """Trees joined onto a parent slot of 'incoming'; own open slot is its dual."""
    if mask.bit_count() == 1:
        j = mask.bit_length()-1
        return (j,) if LABELS[j] == incoming else ()
    result = []
    own = DUAL[incoming]
    for v in VERTICES:
        if own not in v:
            continue
        other = list(v)
        other.remove(own)
        for subsets in parts(mask, len(other)):
            for assignment in sorted(set(permutations(other))):
                for children in product(*(trees(m,t) for m,t in zip(subsets,assignment))):
                    result.append((own, tuple(zip(assignment,children))))
    return tuple(result)


def walk(t):
    if isinstance(t,int):
        return
    yield t
    for _,s in t[1]:
        yield from walk(s)


@cache
def leaves(t):
    if isinstance(t,int):
        return (t,)
    return tuple(i for _,s in t[1] for i in leaves(s))


def kind(t):
    nearh = sum(tuple(sorted((own, *(k for k,_ in kids))))==('F','f','h')
                for own,kids in walk(t))
    h_internal = any(k=='h' and not isinstance(s,int)
                     for _,kids in walk(t) for k,s in kids)
    return ('double' if nearh==2 else
            ('split' if h_internal else 'mixed') if nearh==1 else 'bridge')


def has_chi(t):
    return any(k=='c' for _,kids in walk(t) for k,_ in kids)


def gtree():
    mask = ((1<<len(LABELS))-1) ^ (1<<ROOT)
    graphs = trees(mask,'g')  # open G tensor leg, contracted with external anti-g
    if not (len(graphs)==len(set(graphs))==27 and Counter(map(kind,graphs))==EXPECT):
        raise AssertionError("unexpected graph census")
    for t in graphs:
        if sorted(leaves(t))!=[0,1,2,4,5]:
            raise AssertionError("bad labelled graph")
        if sum(len(kids)-1 for _,kids in walk(t))!=4:
            raise AssertionError("bad valence sum")
    return graphs


def action_rule(types, p, a, m, H, mu):
    """Open vertex tensor: coefficient of i L via square-free multilinear product."""
    shape = tuple(4 if t not in 'ch' else 1 for t in types)
    ans = np.zeros(shape, dtype=complex)
    if np.max(np.abs(sum(p))) > 1e-9:
        raise AssertionError("vertex momentum routing")
    b = mu*a/m
    v = m/a
    if any(t in 'FfGg' for t in types):
        slot = {t:i for i,t in enumerate(types)}
        if 'F' in slot:
            mats = {'w':-1j*a*np.einsum('ij,jkl->ikl', ETA, GAM@PL),
                    'h':(-1j*b*I)[None,:,:],
                    'c':(b*G5)[None,:,:]}
            upper, lower = 'F', 'f'
        else:
            mats = {'w':-1j*a*np.einsum('ij,jkl->ikl', ETA, GAM)}
            upper, lower = 'G', 'g'
        k, = [t for t in types if t not in 'FfGg']
        return np.transpose(mats[k], tuple((k,upper,lower).index(t) for t in types))
    # Explicit derivatives of Fourier-space neutral action; slots of identical
    # fields contribute separately (e.g. 3! terms for hhh, 2! * 2! for wwhh).
    coeff = {
        'hww': ((('h','w','w'), 1j*a*m),),
        'chw': ((('h','c','w'), a), (('c','h','w'),-a)),
        'cch': ((('h','c','c'), -1j*H*H/(2*v)),),
        'hhh': ((('h','h','h'), -1j*H*H/(2*v)),),
        'hhww': ((('h','h','w','w'), 1j*a*a/2),),
        'ccww': ((('c','c','w','w'), 1j*a*a/2),),
        'hhhh': ((('h','h','h','h'), -1j*H*H/(8*v*v)),),
        'cchh': ((('c','c','h','h'), -1j*H*H/(4*v*v)),),
        'cccc': ((('c','c','c','c'), -1j*H*H/(8*v*v)),),
    }[''.join(sorted(types))]
    for monomial,coupling in coeff:
        for order in permutations(range(len(types))):
            if tuple(types[i] for i in order)!=monomial:
                continue
            ws = [i for i,t in enumerate(types) if t=='w']
            for index in np.ndindex(shape):
                factor = 1
                if len(ws)==2:
                    factor = ETA[index[ws[0]],index[ws[1]]]
                elif len(ws)==1:
                    factor = (ETA@p[order[1]])[index[ws[0]]]
                ans[index] += coupling*factor
    return ans


def expanded_sum(graphs, P, external, a,m,H,mu,xi):
    """Recursive raw tensor contractions, no Ward replacements."""
    @cache
    def momentum(t):
        return sum(P[j] for j in leaves(t))
    @cache
    def contract(t):
        own,kids=t
        T=action_rule((own,*(k for k,_ in kids)),
                      (-momentum(t),*(momentum(s) for _,s in kids)), a,m,H,mu)
        terms=[T,list(range(T.ndim))]
        for i,(k,s) in enumerate(kids,1):
            if isinstance(s,int):
                wave=external[s]
            else:
                q=momentum(s)
                q2=dot(q,q)
                if k=='w':
                    prop=-1j*(ETA-(1-xi)*np.outer(q,q)/(q2-xi*m*m))/(q2-m*m)
                elif k=='c':
                    prop=np.array([[1j/(q2-xi*m*m)]])
                elif k=='h':
                    prop=np.array([[1j/(q2-H*H)]])
                else:
                    mf=mu if k in 'fF' else 0.
                    if k in 'fg':
                        prop=1j*(slash(q)+mf*I)/(q2-mf*mf)
                    else:
                        prop=1j*(slash(-q)+mf*I).T/(q2-mf*mf)
                wave=prop@contract(s)
            terms.extend((wave,[i]))
        return np.einsum(*terms,[0])
    return np.array([external[ROOT]@contract(t) for t in graphs])


def unitary_closed(params,data,e):
    """Independent contracted amplitude using only unitary W and scalar Higgs.
    No graph generator, tensor action, Goldstone or common sewing functions."""
    a,m,H,mu,q,K,p,r,_=params
    f,g,u,vb=data
    v=m/a
    b=mu/v
    j=np.array([vb@gg@PL@u for gg in GAM])
    S=lambda l:(slash(l)+mu*I)/(dot(l,l)-mu*mu)
    push=lambda k,J:-1j*(J-k*dot(k,J)/(m*m))/(dot(k,k)-m*m)
    end=-a*e/(dot(K,K)-m*m)
    initial=push(q,-1j*a*j)
    hh=p+r
    split_factor=3*H*H/(v*(dot(hh,hh)-H*H))
    wh=2j*a*m
    bsum = dot(initial,end)*(2j*a*a+wh*split_factor)
    msum = 0j
    for x in (p,r):
        k=q+x
        bsum += wh*dot(push(k,wh*initial),end)
        xspin=np.array([vb@(gg@PL@S(f+x)+S(f-k)@gg@PL)@u for gg in GAM])
        msum += wh*dot(push(k,-1j*a*b*xspin),end)
    W=slash(end)@PL
    ds=0j
    for x,y in ((p,r),(r,p)):
        ds += vb@(W@S(f+x+y)@S(f+x)
                +S(f+x+K)@W@S(f+x)
                +S(f+K+x)@S(f+K)@W)@u
    ds *= -1j*a*b*b
    vs=-1j*a*b*split_factor*(vb@(W@S(f+hh)+S(f+K)@W)@u)
    return dict(bridge=bsum,mixed=msum,double=ds,split=vs)


def accepted(params, P, graphs):
    a,m,H,mu,*_=params
    den=[]
    for t in graphs:
        for _,kids in walk(t):
            for k,s in kids:
                if isinstance(s,int): continue
                q=sum(P[j] for j in leaves(s)); q2=dot(q,q)
                masses=(m*m,*(xi*m*m for xi in GAUGES)) if k=='w' else (
                    tuple(xi*m*m for xi in GAUGES) if k=='c' else
                    ((H*H,) if k=='h' else ((mu*mu,) if k in 'fF' else (0.,))))
                den.extend(abs(q2-ms) for ms in masses)
    return float(min(den))


def exact_identities():
    """Contract the handwritten raw-rule polynomials without kinematic sampling."""
    a,m,M,s,A,B,J,E,X,eQ,xi,u,P,Ph=sp.symbols(
        'a m M s A B J E X eQ xi u P Ph', nonzero=True)
    D=s-m*m
    C=s-xi*m*m
    # A=q², B=q·Q, J=q·e. Contract the W-root rule with E=q.
    qT=4*a*a*m*m*(J-(1-xi)*B*eQ/C)/D-2*a*a*(2*B-A)*eQ/C
    mG=2*a*a*m*m*(2*J-eQ-(1-xi)*(2*B-s)*eQ/C)/D+2*a*a*M*M*eQ/C
    assert sp.simplify((qT-mG+2*a*a*eQ).subs(B,(s+A-M*M)/2))==0
    # Test variation on a GENERIC root E-vector, using E=E·q, X=E·Q.
    T=4*a*a*m*m*(J-(1-xi)*X*eQ/C)/D-2*a*a*(2*X-E)*eQ/C
    coef=2*a*a*m*m*(xi-1)*eQ/(D*C)
    assert sp.simplify(T-T.subs(xi,1)-coef*E)==0
    # Generic root sewing is written in the normalized phase of the derivation.
    Tr,V,Gr=sp.symbols('Tr V Gr')
    source= (Tr-(1-xi)*u*V/C)/D-u*Gr/(m*C)
    assert sp.simplify((source-(Tr-u*V/(m*m))/D).subs(Gr,V/m))==0
    # Mixed source endpoint u=mu P and Q.J_h=mu Ph+P.
    mu=sp.Symbol('mu',nonzero=True)
    b=mu*a/m
    Jr=mu*Ph+P
    right=2*a*m*eQ
    mixed=(b*Tr-b*(1-xi)*Jr*right/C)/D-b*mu*Ph*right/(m*m*C)
    assert sp.simplify(mixed-(b*(Tr-Jr*right/(m*m))/D+b*P*right/(m*m*C)))==0
    assert sp.simplify((mixed-mixed.subs(xi,1)).subs(P,u/mu)-coef*u/(m*m))==0
    return dict(born_ordered_ST_remainder='-2 a^2 e.Q',
                born_root_after_seagull='0; e.(Q_p+Q_r)=e.q',
                born_inner_variation='-2 a^2 u (xi-1) e.Q / (D C)',
                mixed_inner_variation='+2 a^2 u (xi-1) e.Q / (D C)',
                own_far_longitudinal_factor='K.e = 0 identically')

def solve(seed,count):
    identities=exact_identities()
    graphs=gtree()
    kinds=np.array([kind(t) for t in graphs])
    chis=np.array([has_chi(t) for t in graphs])
    rms = {}
    def obs(k,x):
        m=float(np.max(np.abs(x)))
        if not np.isfinite(m):
            raise AssertionError(k+" nonfinite")
        rms[k]=max(rms.get(k,0.),m)
    rng=np.random.default_rng(seed)
    out=[]
    rejects=0
    for _ in range(count):
        for _trial in range(10000):
            found=fixture.get_fixture(rng)
            if found is not None:
                params,e,near,far=found
                a,m,H,mu,q,K,p,r,_=params
                f,F,u,vb=near
                g,G,ug,vbg=far
                P=np.array((f,F,g,G,p,r))
                gap=accepted(params,P,graphs)
                if gap>GAP:
                    break
            rejects+=1
        else:
            raise RuntimeError("fixture retry limit")
        ext=(u,vb,ug,vbg,np.ones(1),np.ones(1))
        un=unitary_closed(params,near,e)
        val=np.array([expanded_sum(graphs,P,ext,a,m,H,mu,xi) for xi in GAUGES])
        swp=np.array([expanded_sum(graphs,P[[0,1,2,3,5,4]],ext,a,m,H,mu,xi) for xi in GAUGES])
        sm=val.sum(axis=1)
        scale=max(1.,np.max(np.abs(val)),np.max(np.abs(sm)),abs(sum(un.values())))
        cm=fixture.cm
        obs('physical_total_nonzero',sm/scale)
        unit_diff=sm-sum(un.values())
        obs('signed_total_vs_independent_unitary_scaled',unit_diff/scale)
        obs('physical_xi_variation_scaled',(sm-sm[2])/scale)
        obs('Bose_total_exchange_scaled',(sm-swp.sum(axis=1))/scale)
        exchange=val[:,kinds=='bridge'].sum(axis=1) + val[:,kinds=='mixed'].sum(axis=1)
        obs('bridge_plus_mixed_vs_unitary_scaled',(exchange-un['bridge']-un['mixed'])/scale)
        for name in ('double','split'):
            obs(name+'_independent_signed_scaled',(val[:,kinds==name].sum(axis=1)-un[name])/scale)
        for name, keep in {'omit_mixed':kinds!='mixed','omit_mixed_chi':~((kinds=='mixed')&chis)}.items():
            ctl=val[:,keep].sum(axis=1)
            obs(name+'_nonzero_scaled',(ctl-ctl[2])/scale)
        obs('external_checks',np.r_[P.sum(axis=0),dot(f,f)-mu*mu,dot(F,F)-mu*mu,dot(g,g),dot(G,G),
                                  dot(K,e),dot(p,p)-H*H,dot(r,r)-H*H,
                                  (slash(f)-mu*I)@u,vb@(slash(F)+mu*I),slash(g)@ug,vbg@slash(G)])
        out.append(dict(parameters=dict(a=a,m=m,M=H,mu=mu),external_momenta=cm(P),
                        external_spinors=cm([u,vb,ug,vbg]),minimum_denominator=gap,
                        comparison_scale=float(scale),
                        unitary_classes={k:cm(v) for k,v in un.items()},
                        results=[dict(xi=xi,signed_graphs=cm(row),signed_total=cm(s),
                                      lhs_minus_unitary=cm(diff),
                                      bose_swapped_total=cm(bose.sum()),
                                      omit_mixed=cm(row[kinds!='mixed'].sum()),
                                      omit_mixed_chi=cm(row[~((kinds=='mixed')&chis)].sum()))
                                 for xi,row,s,diff,bose in zip(GAUGES,val,sm,unit_diff,swp)]))
    for k,v in rms.items():
        if k.endswith('_nonzero'):
            assert v>NEG, (k,v)
        elif k.endswith('_nonzero_scaled'):
            assert v>NEG, (k,v)
        else:
            assert v<EOS, (k,v)
    return dict(seed=seed,accepted=count,rejected=rejects,xi=GAUGES,
                exact_identities=identities, script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                graph_counts=EXPECT, number_goldstone_graphs=int(chis.sum()),
                graph_signatures=[dict(id=i, class_name=kind(t), goldstone=has_chi(t), labelled_tree=t)
                                  for i,t in enumerate(graphs)],
                max_metrics=rms, fixtures=out,
                fixture_source_sha256=PIN)


def main():
    pa=argparse.ArgumentParser(description=__doc__)
    pa.add_argument('--seed',type=int,required=True)
    pa.add_argument('--points',type=int,required=True)
    pa.add_argument('--output',type=Path,required=True)
    a=pa.parse_args()
    if a.points<1: pa.error('--points must be positive')
    result=solve(a.seed,a.points)
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(f"seed={a.seed} accepted={a.points} rejected={result['rejected']}")
    print('graph census', result['graph_counts'], 'Goldstone graphs',result['number_goldstone_graphs'])
    for k,v in result['max_metrics'].items():
        print(f'{k}: {v:.17g}')


if __name__=='__main__':
    main()

---
title: "Boundary data for fixed-root scalar shifts in broken SU(2)"
subtitle: "Theorem, six-point amplitude, and independent contour checks"
date: "25 November 2026"
geometry: margin=24mm
fontsize: 11pt
---

## Abstract

We study a null shift of two equal-mass Higgs momenta that leaves both virtual vector momenta fixed. For a neutral vector bridge sewn to two arbitrary conserved sources in a single-doublet, renormalizable broken \(SU(2)\) theory, the complete reduced tree amplitude is bounded but generically nonvanishing at infinity. We give its boundary as an explicit bilinear in the currents and show that no generic shift direction in this class removes it for every allowed pair of currents. This obstruction persists even when the unphysical \(R_\xi\) poles cancel exactly. We then obtain the large-shift constant for a complete Yukawa-coupled six-point tree, separating the bridge, single emission, double emission from the same fermion line, and stationary scalar splitting. Six moving poles plus the independently derived boundary reconstruct the benchmark amplitude. The result identifies the extra data needed under the fixed-root deformation, not a failure of all-line transverse recursion or a definition of a physical spin density away from a pole.

## Setup and prior work

Metric \(g=\mathrm{diag}(1,-1,-1,-1)\), all momenta incoming, Fourier waves \(e^{-ipx}\). The scalar doublet has a radial Higgs \(h\) of mass \(M\), vacuum expectation value \(v\), and a neutral \(W^3\) of mass \(m=av\), with \(a=g_2/2\). Write \(q,K\) for the two end momenta and \(p,r\) for the two Higgs momenta, \(q+K+p+r=0\), \(p^2=r^2=M^2\). Internal exchanges have propagators \(1/(Q^2-m^2)\) and \(1/(Q^2-\xi m^2)\); no widths are used. A dot product is bilinear, with no complex conjugation. The Higgs deformation is
\[
p(z)=p+zn,\qquad r(z)=r-zn,\qquad
n^2=n\cdot p=n\cdot r=0,\qquad d=2q\cdot n\ne0.
\]
The test keeps \(q,K,p+r\) fixed and excludes coincident moving poles and stationary singularities. In the Higgs-pair rest frame the two projective null choices are the complex linear combinations of the two spatial directions transverse to the Higgs axis. Our result applies to both choices when the nondegeneracy conditions hold.

Ke, Liu and Mahbub derive an off-shell continuation of the on-shell transverse/longitudinal projectors and discuss the large-\(z\) conditions for the **all-line transverse (ALT)** shift [arXiv:2611.13096v1, §§3–4 and App. B](https://arxiv.org/html/2611.13096v1#S3). Their constructibility statements for ALT do not imply that this two-Higgs shift has a vanishing boundary. Our earlier [fixed-root check](https://github.com/lm846802297-bit/off-shell-vector-polarization/tree/6d92e34be9694920dc0c64fd0a1d18747ccb6622/research-baselines/su2-fivepoint-fixedq-boundary) tested the off-shell *contracted* weak current against the Goldstone current. Here we treat the physical contraction between two sources and then the full signed six-point amplitude. We take the on-shell pole factorization as given; the added result is the physical boundary in this deformation class.

## Two-source boundary theorem

Assume specified sources \(j,e\) with \(q\cdot j=K\cdot e=0\), which are fixed under the shift and couple to the neutral bridge without separately emitting \(h\). This includes two distinct zero-Yukawa fermion currents in the stated neutral channel. The fixed root propagators and overall phase are removed from the reduced amplitude \(\mathcal B\). Let \(Q_\pm=q+p,q+r\) at \(z=0\), \(D_\pm=Q_\pm^2-m^2\), and \(s_{hh}=(q+K)^2\). The sums of the two derivative Goldstone vertices and the gauge-dependent part of the \(R_\xi\) vector propagator reduce the ordered numerator to the unitary expression, since
\[
\frac{1+m^2(1-\xi)/(Q^2-m^2)}{Q^2-\xi m^2}
=\frac{1}{Q^2-m^2}.
\]
Adding the \(WWhh\) vertex and the stationary \(h^*\to hh\) split then gives
\[
\mathcal B(z)=
C_0\, j\cdot e+4a^2\sum_{\sigma=\pm1}
\frac{m^2 j\cdot e-[j\cdot(Q_\sigma+\sigma zn)]
[e\cdot(Q_\sigma+\sigma zn)]}{D_\sigma+\sigma dz},
\quad C_0=2a^2\left(1+\frac{3M^2}{s_{hh}-M^2}\right).
\]
The products quadratic in \(z\) are the same in the two orderings. Dividing by their denominators, the terms linear in \(z\) cancel; the remaining constant is
\[
B_\infty=C_0\,j\cdot e
-4a^2\left[
  \frac{[j\cdot(q-K)](e\cdot n)+(j\cdot n)[e\cdot(q-K)]}{d}
  -\frac{(j\cdot n)(e\cdot n)(D_++D_-)}{d^2}
\right].
\]
Here \(D_++D_-=2(M^2-m^2-q\cdot K)\). The expression is unchanged by rescaling \(n\), but can change between the two projective directions. Given the simple finite residues \(\rho_\sigma\) in the \(z\) plane, the exact recursive identity for this reduced amplitude is
\[
\mathcal B(0)=B_\infty-\sum_\sigma\frac{\rho_\sigma}{z_\sigma}.
\]
The scalar-splitting term in \(B_\infty\) has an \(s_{hh}=M^2\) pole in the full amplitude. Calling all of \(B_\infty\) a local contact interaction would be incorrect: it is regular in \(z\), not necessarily in the unshifted invariants.

**Obstruction valid for a class of sources.** At any admissible \(n\),
the Gram determinant of \(q,K,n\) is \(- (q\cdot n)^2 s_{hh}\).
If it is nonzero, the one-dimensional common orthogonal complement is
nondegenerate. Choose \(j=e=t\) in that complement, \(t^2\ne0\).
All ordered boundary terms vanish and \(B_\infty=C_0 t^2\ne0\) for
generic \(s_{hh}\). Thus no admissible direction at generic kinematics
eliminates the boundary as a bilinear on the whole conserved-current
space. This is not a claim that every individual helicity choice has a
nonzero boundary. For timelike roots, varying the directions of a
massless fermion pair spans the transverse current space; the witness
can be understood as a linear combination of currents rather than a
single real helicity state.

## Complete six-point benchmark with nonzero Yukawa coupling

We reuse the [locked 27-graph amplitude](https://github.com/lm846802297-bit/off-shell-vector-polarization/tree/ba7dfb747da24809d1c1d9250bd9f8f6ce3677bf/research-baselines/su2-sixpoint-connected-tree) for distinguishable fermion pairs. The \(f\) line has \(f_L\) in a doublet and \(f_R\) a singlet with Yukawa mass \(\mu\); the far line is massless and vectorlike. The near chiral current \(j\) need not be conserved, while \(K\cdot e=0\). Strip no factors in this section. In unitary gauge define
\[
L=-a\,\frac{j-q(q\cdot j)/m^2}{q^2-m^2},\qquad
R=-a\,\frac{e}{K^2-m^2},\qquad
W=\not R P_L,\quad b=\frac{\mu a}{m}.
\]
In the signed amplitude convention of the benchmark, the bridge
boundary is \(iB_\infty[L,R]\): the fixed root factors are now included
in \(L,R\), and the same polynomial division needs no additional
conservation condition in unitary gauge.

For the remaining classes, the hard fermion propagator obeys
\[
S(P+wn)=\frac{\not n}{\ell}+\frac{1}{w}
  \left(\frac{\not P+\mu}{\ell}
       -\frac{(P^2-\mu^2)\not n}{\ell^2}\right)
  +O(w^{-2}),\qquad \ell=2P\cdot n.
\]
Here \(w\) is \(+z\) or \(-z\), and all hard-line slopes are taken nonzero.
For a Higgs assignment \(x=p\) or \(r\), set \(k=q+x\) and define \(J_0,J_1\) as the first two Laurent coefficients of
\[
J^\alpha(w)=\bar v_f[
      \gamma^\alpha P_L S(f+x+wn)
     +S(f-k-wn)\gamma^\alpha P_L]u_f .
\]
The coefficient \(n\cdot J_0\) vanishes by \(\not nP_L\not n=0\); the leading longitudinal numerator in the internal unitary propagator therefore does not produce a growing term in this class. Writing \(g_h=2iam\), its constant is
\[
B_{\rm mixed}=\sum_{x=p,r}
 \frac{g_h a b\,(R\cdot n)}{2m^2 k\cdot n}\,
 (k\cdot J_0+n\cdot J_1).
\]
For the double line emission, write \(U=\not n/(2f\cdot n)\) and
\(V=\not n/[2(f+K)\cdot n]\). The two label orderings give the same
leading constant, since \(\not n/\ell\) is unchanged under \(n\to-n\):
\[
B_{\rm double}=-2iab^2\bar v_f\big[
 W S(f+p+r)U+V W U+V S(f+K)W\big]u_f .
\]
The line-emitted scalar split is stationary and its value is its boundary:
\[
B_{\rm split}=-iab\frac{3M^2}{v(s_{hh}-M^2)}
 \bar v_f\big[W S(f+p+r)+S(f+K)W\big]u_f .
\]
These four classes exhaust the locked benchmark. In all classes the
large-\(z\) amplitude is bounded; the formula makes no assertion of
separate gauge invariance of the classes at finite \(z\).

## Numerical evidence and what remains open

The two-source formula was checked by separate finite \(R_\xi\) contractions,
small contours at the physical and would-be gauge poles, and two outer
circles. In 72 evaluations (including explicit nonzero witnesses)
the largest scaled error in the boundary was \(1.7\times10^{-14}\),
the largest physical-residue mismatch \(6.2\times10^{-14}\), and the
largest unphysical residue \(7.5\times10^{-15}\).

For the full Yukawa calculation, all 20 pinned fixtures and both null
directions independently gave six distinct moving poles, with nonzero
\(B_\infty\) in all 40 cases. Extracting the residues by small circles
and restoring the analytic boundary gave a maximum scaled reconstruction
error of \(1.1\times10^{-13}\). On four fixtures, both shifts and three
\(\xi\) values, the separate Fourier-action graph contraction confirmed
the outer-contour boundary within \(7.5\times10^{-11}\). This checks
the signs of the four calculated pieces without using gauge identities
to simplify the test graph sum. Errors are numerical residuals, not
theoretical uncertainties.

This result does not cover arbitrary broken gauge groups, non-Abelian
charged channels, all choices of shifted legs, ALT, loop amplitudes,
widths or complex physical poles. It does not define a gauge-independent
helicity fraction for a generic off-shell current. The natural next
test for a more general article is to determine how the witness and
hard-line cancellation change when the neutral commuting bridge is
replaced by a charged bridge with genuine gauge self-interactions.
Before any submission, compare against the closest related versions in
the literature and have a human author reproduce the amplitude identities.

## Reproduction

The accompanying Python files provide the closed formulas and the
independent contour checks. The second script imports the previously
archived six-point fixture/evaluator, identified by commit
\texttt{ba7dfb747da24809d1c1d9250bd9f8f6ce3677bf}. The archive includes
the exact hash of each parent file, both raw outputs, and the commands.

# When do moving poles suffice? The non-Abelian \(W^+W^-W^3hh\) tree

## Result and conventions

For three physical equal-mass vectors in a renormalizable single-doublet
\(SU(2)\) theory, the complex two-Higgs shift gives a bounded, generically
nonzero tree amplitude. The exact boundary separates into (i) the seagull
and stationary Higgs exchange, (ii) ordered emissions on a single vector
arm, and (iii) emissions from different arms of a true cubic gauge vertex.
The first two use the same fixed-root building blocks as the
[neutral result](../su2-fixed-root-boundary/OffShell_Polarization_Boundary.md).
The third cannot be obtained by sewing just one neutral two-source boundary.
We give a computable zero-boundary criterion and a nonzero stationary-pole
witness. This is a result for the specified shift and field content, not
a statement about all electroweak amplitudes.

All momenta are incoming, metric \(+---\); \(k_i^2=m^2\),
\(e_i\cdot k_i=0\), \(p^2=r^2=M^2\), and \(k_1+k_2+k_3+p+r=0\).
Here \(m\) is the common vector mass, \(M=m_h\), \(a=g_2/2\), and \(v=m/a\).
Work in unitary gauge and factor out a common cubic-group factor and
overall phase. In the real adjoint basis one may use legs \(1,2,3\);
changing to \(+,-,3\) only changes the common color phase, because each
surviving tree contains exactly one antisymmetric cubic group tensor,
and scalar vertices and propagators are color diagonal. The cubic
kinematic vertex is
\[
V(k_1,k_2,k_3;e_1,e_2,e_3)
=\sum_{\rm cyc}(e_i\cdot e_j)(k_i-k_j)\cdot e_s.
\]
Use \(U(Q)e=[e-Q(Q\cdot e)/m^2]/D(Q)\), \(D(Q)=Q^2-m^2\),
\(g_h=2am\), \(g_{hh}=2a^2\).
Let \(H=p+r\) and shift \(p\mapsto p+zn,\ r\mapsto r-zn\), with
\(n^2=n\cdot p=n\cdot r=0\). Require \(\ell_i=2k_i\cdot n\ne0\) and
nondegenerate propagators; complex contractions are bilinear with no
conjugation. The space of projective null directions has two choices
at generic equal-mass \(p,r\) kinematics.

## Complete graph count and full amplitude

A five-leaf tree has either three cubic vertices or one cubic and one
quartic. The only vertices of unitary-gauge single-doublet \(SU(2)\)
are \(WWW\), \(WWWW\), \(WWh\), \(WWhh\), \(hhh\), \(hhhh\).
An odd number of external vectors requires a cubic \(WWW\).
There are **18 labelled trees**:

- 6 ordered emissions of both Higgs particles on a single arm;
- 6 placements on different arms;
- 3 \(WWhh\) seagulls;
- 3 \(h^*\to hh\) splits.

The quartic gauge vertex cannot attach to two Higgs particles through
one allowed cubic vertex. A separate partition recursion with all these
vertex types and fixed external adjoint colors finds the same count.

For each arm \(i\), define \(K_i=k_i+H\) and
\[
 C_0=2a^2\left(1+\frac{3M^2}{H^2-M^2}\right),\qquad
 E_i(p,r)=U(K_i)\{C_0e_i+g_h^2[U(k_i+p)+U(k_i+r)]e_i\}.
\]
A term with both Higgs particles on arm \(i\) is \(V\) with the \(i\)-th
arguments replaced by \(K_i,E_i\). For distinct arms \(i,j\), replace
those two arguments by \(k_i+x,g_hU(k_i+x)e_i\) and
\(k_j+y,g_hU(k_j+y)e_j\), summing over cyclic arm pairs and \((x,y)=(p,r),(r,p)\).
This is the complete signed amplitude before taking any large-\(z\) limit.

## Closed boundary, and a practical cancellation test

Set \(R_i=(n\cdot e_i)/\ell_i\) and
\[
 F_i=\frac{H\cdot e_i-2R_i(M^2+k_i\cdot H)}{\ell_i},\qquad
 S_i=-\frac{R_i(2k_i+H)+nF_i}{m^2}.
\]
For the same-arm terms, expand \(U(k_i+x+wn)e_i\) for \(w=\pm z\).
Its linear part is \(-w\,nR_i/m^2\), which cancels between the two
orders. The sum of constants is \(S_i\). Thus the same-arm plus fixed
part \(B_{\rm same}\) is obtained by inserting
\[
 B_i=U(K_i)(C_0e_i+g_h^2 S_i)
\]
in the \(i\)-th arm of \(V\), again using momentum \(K_i\).
The other two physical equal-mass vectors yield a conserved \(WWW\)
current at that arm, as follows directly from contracting the cubic
vertex with its root momentum. This is why the neutral fixed-root
two-source construction can be sewn into this part.

For different arms, let \(L_i=k_i+x+wn,\ L_j=k_j+y-wn\),
\(L_i+L_j+k_s=0\), \(b_i=L_i\cdot e_i,\ b_j=L_j\cdot e_j\).
The two contractions with \(L_i\) and \(L_j\), using \(k_s^2=m^2\) and
\(k_s\cdot e_s=0\), yield the exact identity
\[
\begin{split}
 V(L_i,L_j,k_s;U(L_i)e_i,U(L_j)e_j,e_s)
 &= \frac{V(L_i,L_j,k_s;e_i,e_j,e_s)}{D_iD_j}\\
 &+\frac{1}{m^2}\left[\frac{b_i}{D_i}e_j\cdot e_s
   -\frac{b_j}{D_j}e_i\cdot e_s
   +\frac{b_ib_j}{D_iD_j}L_i\cdot e_s\right].
\end{split}
\]
Here \(D_{i,j}=D(L_{i,j})\). The first term decays as \(1/w\);
the last term seems to grow as \(w\), but its growth cancels when
the labelled Higgs assignments are paired. The two single-ratio constants proportional to \(R_i(e_j\cdot e_s)\) and \(R_j(e_i\cdot e_s)\) cancel in the cyclic sum. The surviving cross-arm boundary is
\[
 B_{\rm cross}=\frac{g_h^2}{m^2}\sum_{\rm cyc}
 \left[ R_iR_j(k_i-k_j)\cdot e_s+
 (n\cdot e_s)(R_jF_i-R_iF_j)\right].
\]
There is no growing part of the sum. The boundary is exactly
\(B_\infty=B_{\rm same}+B_{\rm cross}\). Given \((k_i,e_i,p,r,n)\),
these formulas use the unshifted vectors and scalar invariants and no
large-\(z\) fit. The moving physical poles alone recover \(A(0)\) **if
and only if** this sum vanishes at the point. In particular:

- If at least two physical polarizations obey \(n\cdot e_i=0\), then
  \(B_{\rm cross}=0\) identically. One such polarization does not suffice.
  The stationary and same-arm terms need not vanish.
- The *full* vanishing condition is a linear kernel in any one vector
  polarization when the other two are fixed. The explicit controls
  include nonzero amplitudes lying in this kernel, rather than setting
  the entire process to zero. The kernel is over complex physical
  polarizations; it is not in general a statement about a standard
  helicity eigenstate.

The boundary contains a physical stationary pole in \(H^2\), so the full
boundary cannot vanish identically for either regular branch of the
shift. Its residue at \(H^2=M^2\) is
\((3M^2/v) A_4(W^+W^-W^3h)\). The other terms are generically regular
at that pole and cannot cancel it. For a concrete witness, choose
opposite incoming three-momenta along \(z\), outgoing vector and
scalar along \(x\), the first two polarizations along \(y\) and the
third along \(z\). Then \(A_4=2g_h|\mathbf k_1|\sum_i 1/D(K_i)\);
this is not identically zero, as seen also by its large-energy
coefficient \(-3/(4 E^2)\) in the denominators' sum. The physical
five-point \(hh\) region has \(H^2\ge4M^2\); the pole argument uses
analytic continuation, not a claimed real \(h\to hh\) decay.

## Independent checks

The recorded 40 deformations (20 real \(2\to3\) phase-space points,
two null directions) have six distinct moving vector poles. For each
pole the three analytically dual spin states sew an \(WWh\) vertex to
the complete four-point \(WWWh\) tree, and agree with a small-contour
extraction to \(8.9\times10^{-14}\) scaled. The independent rooted
tree evaluator uses explicit color permutations of the derivative
slots, rather than the grouped formula; agreement is
\(2.6\times10^{-15}\). A contour outside all moving poles checks
the constant to \(2.2\times10^{-13}\), and restoring it reconstructs
the unshifted result to \(7.9\times10^{-15}\). The output also includes
positive-power contour checks, two-polarization cancellation controls
that retain a nonzero full boundary, and full-boundary cancellation
controls that retain a nonzero amplitude. No finite-\(\xi\) tensor
cross-check is claimed for the charged process.

## What is and is not new relative to known results

- [Ke–Liu–Mahbub, arXiv:2611.13096v1 §§3–4, Appendix B](https://arxiv.org/html/2611.13096v1#S3)
  derive a polarization continuation from on-shell factorization and use
  ALT for their constructibility result. Our fixed-root scalar shift
  leaves some propagator invariants constant and gives *linear* moving
  denominators; it does not obey the large-\(z\) assumptions of their
  ALT bound. This is not an ALT counterexample or a different physical
  massive spin sum.
- [Ema et al., arXiv:2403.15538, §2.2](https://arxiv.org/html/2403.15538#S2.SS2)
  explicitly distinguish their ALT from null shifts whose propagators
  scale differently, and discuss spin-dependent constructibility.
  Shift dependence itself is not a novelty here. We supply explicit
  gauge/Higgs boundary coefficients and polarization cancellation
  conditions for one different five-point deformation.
- [Franken–Schwinn, arXiv:1910.13407, §5](https://arxiv.org/abs/1910.13407)
  already analyze massive tree constructibility and when more lines
  suffice. Their existence and scaling statements do not identify
  the fixed-vector two-Higgs boundary written above.
- [Durieux, Kitahara, Shadmi and Weiss, JHEP 12 (2019) 058, §5.2](https://doi.org/10.1007/JHEP12%282019%29058)
  explicitly discuss adverse large-\(z\) behavior and reconstruct the
  four-point \(W\overline W hh\) amplitude from all physical invariants
  plus a nonfactorizable term determined by the ultraviolet behavior.
  This is close precedent: neither the existence of a Higgs contact
  nor failure of a two-line shift at four points is our novelty.
  The present result adds the two moving legs of the cubic gauge
  vertex, their cross-arm term and its polarization selection rule,
  and distinguishes poles stationary in \(z\) from interactions local
  in all invariants at five points.
- [Jin–Feng, arXiv:1507.00463](https://arxiv.org/abs/1507.00463)
  and [Huang–Jin–Li, arXiv:2210.07025](https://arxiv.org/abs/2210.07025)
  supply systematic boundary-operator methods. We do not claim to
  invent boundary operators or the rational contour identity. Our
  scoped contribution is the explicit decomposition and the algebraic
  two-polarization selection rule for a genuine broken-gauge five-point
  process, paired with the neutral two-source and Yukawa six-point
  checks.

The case is no longer merely one more example of a nonzero boundary:
it resolves what must be added, what vanishes as a tensor contraction,
and why even the simpler common contribution can survive. It is still
a narrow result. We do not claim general-multiplicity, nondegenerate
electroweak, complex-pole or loop theorems, or guaranteed publication.

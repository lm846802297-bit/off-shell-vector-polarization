# From the known four-point amplitude to the five-point selection rule

## 1. Four-point comparison at identical shift

[Balkin, Durieux, Kitahara, Shadmi and Weiss, JHEP 03 (2022) 129, §3.2.2, Eq. (3.37)](https://doi.org/10.1007/JHEP03(2022)129)
give the full four-point amplitude with two vector-exchange channels,
the Higgs pole and the local contact term. Use the same four diagrams
in a Lorentz-vector representation, dropping a common overall amplitude
phase. The comparison below uses normalized vector states directly; it
does not rely on identifying little-group component labels with a
particular vector helicity convention.

Use the conventions of [the parent derivation](CHARGED_BOUNDARY.md).
Let `k,k',p,r` be incoming, with vector polarizations `u,w`; let
\(L_x=k+x\), \(D_x=L_x^2-m^2\), \(d=2k\cdot n\), \(H=p+r\).
For the ordered exchange of the vector carrying momentum \(L_x\),
\[
 A_x(z)=g_h^2\,
 \frac{u\cdot w-(X_x+\sigma_x z N)(Y_x+\sigma_x z P)/m^2}
      {D_x+\sigma_x zd},
 \quad (X_x,Y_x)=(u\cdot L_x,w\cdot L_x),
 \quad (N,P)=(u\cdot n,w\cdot n),
\]
where \(\sigma_p=+1,\sigma_r=-1\). Long division gives
\[
 \frac{A_x(z)}{g_h^2}=
 -\frac{\sigma_x NP}{m^2d}z+
 \frac{1}{m^2}\left[\frac{NP D_x}{d^2}
                      -\frac{N Y_x+P X_x}{d}\right]+O(z^{-1}).
\]
Thus it is the *sum* of the two ordered exchanges, not either one,
that is bounded at infinity. The two terms stationary in \(z\) have
coefficients \(g_{hh}\) and \(g_h\lambda_{hhh}/(H^2-M^2)\),
respectively, multiplying \(u\cdot w\). Here
\(g_{hh}=2a^2\), \(g_h=2am\), \(\lambda_{hhh}=3M^2/v\).
The latter is a physical pole in \(H^2\), not a local contact.
With \(u\cdot k=0\), the ordered constant contracted with \(w\)
is precisely \(g_h^2 w\cdot S\), where \(S\) is the vector in the
parent derivation. This tensor formula also applies with a conserved
source replacing \(w\); its validity for that source follows from
the displayed propagator numerator, not from uniqueness of off-shell
continuation of the on-shell amplitude.

For five points, denote by \(J_i\) the cubic current made from the
other two physical vector legs. Direct contraction of the cubic vertex
gives \(K_i\cdot J_i=0\) for equal masses. Then
\(J_i\cdot U(K_i)X=J_i\cdot X/D(K_i)\).
All 12 same-arm and stationary diagrams therefore arise by inserting
the four-point *tensor* with \(w=J_i/D(K_i)\).
The six remaining diagrams insert one Higgs on each of two different
arms. The cubic vertex itself then has two shifted vector momenta.
This is not a single insertion of the four-point two-Higgs tensor.
It is also **not a new independent five-point coupling**: the cross
contribution uses two ordinary \(WWh\) vertices and one \(WWW\)
vertex. Here “cross contribution” denotes the topology-defined
decomposition in the stated unitary-gauge representation; we do not
assert gauge invariance of each unsummed class in an arbitrary gauge.

## 2. Exact rule when one polarization is orthogonal to the shift

Suppose \(n\cdot e_s=0\) and \((i,j,s)\) is cyclic.
With \(\ell_i=2k_i\cdot n\), \(R_i=(n\cdot e_i)/\ell_i\), we have
\(R_s=0\), \(F_s=(H\cdot e_s)/\ell_s\).
The two other cyclic terms proportional to \(F_s\) combine with
the \((i,j,s)\) term to give
\[
 B_{\rm cross}=\frac{g_h^2R_iR_j}{m^2}\,
 e_s\cdot\left[k_i-k_j+\frac{\ell_j-\ell_i}{\ell_s}H\right]
 =-\frac{2g_h^2R_iR_j}{m^2\ell_s}\,e_s\cdot T_s,
 \qquad T_s=\ell_j k_i-\ell_i k_j.
\]
The last equality uses \(H=-k_i-k_j-k_s\),
\(k_s\cdot e_s=0\), and \(\ell_i+\ell_j+\ell_s=0\).
For nonzero coupling and \(R_iR_j\ne0\), cancellation is **equivalent**
to \(e_s\cdot T_s=0\). Two orthogonal vectors remain sufficient but
are not necessary. The condition is on the cross term only.

For clarity on its scope, put
\(P_s=T_s-n(T_s\cdot k_s)/(n\cdot k_s)\).
Since \(T_s\cdot n=0\), \(P_s\) belongs to the two-dimensional
complex plane \(k_s^\perp\cap n^\perp\), with
\(e_s\cdot P_s=e_s\cdot T_s\). The plane is nondegenerate because the
Gram determinant of \((k_s,n)\) is \(-(k_s\cdot n)^2\ne0\).
For \(P_s\ne0\) the criterion picks a one-dimensional kernel in the
plane, and for \(P_s=0\) the whole plane works. This is a covariant
criterion for complex polarization combinations, not generically a
statement about one standard helicity or one real linear polarization.
Generic non-null kernel vectors can be normalized. The full
criterion for using only the six moving residues still requires
\(B_{\rm cross}+B_{\rm same}=0\).

## 3. Exact regression tests and manuscript claim

Run `python same_shift_check.py` (SymPy; all comparisons exact).
The four-point fixture has a stationary coefficient \(-59/91\) and
two ordered constants \(41/64\) each, giving \(1843/2912\) in total.
The same-arm sewing is checked separately on all three arms.
For the five-point fixture with unit vector mass and gauge coupling,
the real polarizations \(e_1=e_2=(0,0,0,1)\),
\(e_3=(0,1,0,0)\) give \(B_{\rm cross}=1/4\) after the common
cubic color/phase factor is removed. Only \(e_3\) is shift orthogonal.
Another physical polarization in the same plane for leg 3 gives
exactly zero without making either of the other legs shift orthogonal.
A two-orthogonality case also gives zero. Direct rational expansion
of the original six diagrams reproduces all three outcomes for both
regular null directions, without using the compact boundary formula.
These are exact fixtures and do not replace the analytic argument.

A defensible manuscript claim is the compact five-point cross term
and its necessary-and-sufficient criterion on the one-orthogonal
polarization subspace for the specified theory and shift. The old
four-point amplitude and its contact term are checks, not claims of
novelty. The full boundary is not a local operator inferred from
five external legs, and the analysis makes no claim against
all-line transverse constructibility.

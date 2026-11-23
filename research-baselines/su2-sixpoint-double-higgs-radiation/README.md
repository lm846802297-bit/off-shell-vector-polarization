# Same-line double Higgs radiation and virtual h splitting, with a conserved far source

**Scope and baseline.** The [single-Higgs radiation checkpoint](https://github.com/lm846802297-bit/off-shell-vector-polarization/tree/ae185d2ac82c3d54edef574aeb6711e921abe7f0/research-baselines/su2-sixpoint-single-higgs-radiation) established cancellation between the earlier paired exchange subset and diagrams with one Higgs emitted from the massive line and one from the bridge. Here we check two other specific classes: two direct Higgs insertions on the massive line, and one virtual Higgs from that line that splits to two. These are amplitude computations, not a complete diagram census or a claim of full six-point gauge invariance. The distinction between an off-shell Lagrangian identity and on-shell constructibility is as in Ke–Liu–Mahbub [arXiv:2611.13096v1, §4.2](https://arxiv.org/html/2611.13096v1#S4.SS2).

The model and phases are unchanged: fixed diagonal lower T³ component of the massive left doublet and a right singlet, real Yukawa, distinguishable vectorlike massless current at the far end, metric (+---) and all incoming momenta. Let a=g/2, m=a v, b=μ/v, M=m_h, where μ is the massive fermion mass. The bare vertices are −ia γ^α P_L, −ib, and −i3M²/v for W–f–f, h–f–f, and hhh respectively, under Fourier convention exp(−ip·x). The far current e carries momentum K and satisfies K·e=0. It has no Goldstone–fermion vertex. Common far coupling and isospin factors are stripped. At the near line, momenta are absorbed from the u(p_f) endpoint toward vbar(p_b), so q=p_f+p_b and p_f+p_b+K+p+r=0.

## The six direct chains and two splitting chains

Write S(l)=(slash(l)+μ)/(l²−μ²), with bare propagator iS(l). For fixed Higgs order x then y along the fermion flow, the three direct near-line matrices, with b² and the common −ia removed, are

```text
J_xy,1^alpha = vbar S(p_f+K+x) S(p_f+K) γ^alpha P_L u
J_xy,2^alpha = vbar S(p_f+x+K) γ^alpha P_L S(p_f+x) u
J_xy,3^alpha = vbar γ^alpha P_L S(p_f+x+y) S(p_f+x) u .
```

Sum these for both (x,y)=(p,r) and (r,p); this gives six orderings with no extra identical-particle factor for the labelled amplitude. The two splitting orders replace x,y by H=p+r in a single Yukawa insertion. Their two near matrices, stripped of b and −ia, are vbar S(p_f+K) γ^alpha P_L u and vbar γ^alpha P_L S(p_f+H)u. Multiply their sum by the signed splitting factor B=3M²/[v(H²−M²)], from the scalar propagator i/(H²−M²) and bare hhh vertex. The two normalized source coefficients are b² for the six direct orders and bB for the two splitting orders.

Within these two classes the single bridge toward the far pair is a W carrying K. A Goldstone has zero coupling at that far end under the stated assumptions, even though a massive source could emit it at the near end. The tree table here counts only these classes.

## A nontrivial near-current test and far projection

When K enters the W vertex on the near line and R=L+K for the oriented fermion momenta after and before that vertex,

```text
slash(K) P_L = S^{-1}(R) P_L − P_R S^{-1}(L) − μ γ5 .
```

Between the on-shell spinors, the inverse propagators telescope through the three W positions. Denote by J_xy the sum of the three displayed matrices and by P_xy the same chains with γ^alpha P_L replaced by γ5. At fixed x then y,

```text
K·J_xy = −μ P_xy − vbar[S(p_f+x+K)γ5 + γ5 S(p_f+x)]u .
```

The leftover terms come from P_L−P_R=−γ5 next to the Yukawa vertices. Using p_f+x+K=−p_b−y gives the equivalent two-contact expression after summing the Higgs orders. For the two virtual-H orders, K·J_H=−μP_H−vbar γ5 u. This makes clear why near-current conservation must not be imposed simply to remove gauge dependence.

Instead only the far current is conserved. For either normalized class J, direct contraction of the W propagator gives

```text
A_ξ = [J·e − (1−ξ)(K·J)(K·e)/(K²−ξm²)]/(K²−m²)
    = J·e/(K²−m²).
```

Neither remaining class has a new ξ variation for the physical far source, though both have generally nonzero finite amplitudes. For the intentionally wrong far source e'=e+αK, the change relative to ξ=1 is (ξ−1)(K·J)(K·e')/[(K²−m²)(K²−ξm²)], generally nonzero in both classes. This is a numerical diagnostic, not an extension of the physical model.

## Independent checks and replay

The numerical route independently walks exact ordered bare matrix chains and uses the explicit R_ξ propagator. It generates the h vertex by a derivative of the Yukawa term and the hhh vertex by Fourier coefficient extraction from the scalar action. It compares every raw amplitude with the separately grouped formula, checks both endpoint identities and Bose exchange of the exact ordered blocks, and checks physical and wrong-source gauge variations with explicit on-shell spinor solutions. The external-wavefunction fixture and scalar action helper are reused (not independent physics evidence for the preceding class) from [the exact earlier code](https://github.com/lm846802297-bit/off-shell-vector-polarization/blob/ae185d2ac82c3d54edef574aeb6711e921abe7f0/research-baselines/su2-sixpoint-single-higgs-radiation/check_higgs_emission.py), blob 25d04b734a686f16cbc090cb89e67b024bfa4373, SHA-256 d05358a11f698dea04e34cf60b5017c42de3d6e6f131df8012f8855a941e62d2. The raw ordered chains, hhh splitting contraction, and endpoint tests are new. Across the 20 accepted fixtures at ξ=0.43, 1, 2.65, the largest physical relative gauge variation was 1.62e-16, per-order comparison 7.61e-17, and wrong-source comparison 1.95e-16. Both wrong-source gauge variations were nonzero, reaching 0.198 and 0.122 in the archived scaling. The finite physical amplitudes of both classes are likewise nonzero; these residuals are numerical checks, not theoretical uncertainties.

```sh
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 python check_double_higgs.py --seed 132023 --points 8 --output results-primary.json
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 python check_double_higgs.py --seed 629431 --points 12 --output results-holdout.json
sha256sum -c SHA256SUMS
```

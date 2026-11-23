# Complete connected neutral SU(2)–Higgs six-point tree

The fixed-component tree $f\bar f,g\bar g,hh$ has **27 labelled graphs**: 11 with both Higgs bosons emitted by the bosonic bridge, 8 with one Higgs off the massive fermion line, 6 with both Higgs bosons emitted directly from that line, and 2 with a virtual Higgs emitted there and split to $hh$. Their full signed sum is independent of $\xi$ and equals an independently assembled unitary-gauge amplitude. Neither the 11 nor 8 graph subset is independently gauge invariant. This is one on-shell-external scattering amplitude, not a proof for arbitrary off-shell currents, isospin choices, recursion, loops, widths, or complex poles.

## Model and exact census

Gauge group SU(2), one scalar doublet $\Phi=(\phi^+,(v+h+i\chi)/\sqrt{2})$, (a=g_{SU(2)}/2, m=av, M=m_h). The near (f_L) is the lower component of a left doublet, (f_R) is a singlet, and $\mu=yv/\sqrt2$ comes from $-y\bar Q_L\Phi f_R+\mathrm{h.c.}$. The far (g_L,g_R) are the lower components of a separate massless vectorlike doublet, with no Yukawa or mixing with the near species. The lower diagonal (T^3=-1/2) convention is $D=\partial+i a W$ on the lower scalar, with (W=W^3). Potential $V=M^2(2vh+h^2+\chi^2)^2/(8v^2)$ in the neutral restriction, expanded at a stationary tree vacuum. Any anomaly-cancelling spectators do not join these external fermion lines at tree level.

Although no local SU(2) generator leaves the vev invariant, the above interactions preserve the **global** charge (T^3+X), where (X=+1/2) on the scalar and both matter doublets, and (X=0) on the near right singlet. The vev and lower components are neutral; $W^\pm,\phi^\pm$ and the upper matter components carry nonzero charge. The linear $R_\xi$ gauge fixing $\mathcal L_{gf}=-(\partial\cdot W-\xi m\chi)^2/(2\xi)$ in the stated neutral convention cancels quadratic $W$–$\chi$ mixing and extends covariantly under this global charge to the charged components. The tree vacuum has no tadpole. Each Faddeev–Popov vertex carries one ghost and one antighost; without external ghosts such internal lines would form a directed cycle and cannot enter a tree.

Each internal edge of a connected tree is a cut. A bosonic cut cannot split either fermion pair: the cut carries no near or far fermion number, separately conserved. Its external subset is a union of whole pairs and Higgs leaves, hence has no global charge. This excludes a charged boson **everywhere**, not just on the bridge. A matter cut must carry its chain's species and external component charge; the upper component cannot occur. Consequently the neutral vertices are exhaustive. For fermions they are $f\bar f W,f\bar f\chi,f\bar f h,g\bar g W$. For bosons they follow from

```text
L3 = a W·(h ∂χ − χ ∂h) + a m W² h − M² h(h²+χ²)/(2v),
L4 = (a²/2) W²(h²+χ²) − M²(h²+χ²)²/(8v²).
```

There is no $Whh$ or neutral $WWW$. Every neutral bosonic vertex is even under $(W,\chi)\mapsto(-W,-\chi)$. After contracting each fermion chain to one node there is one connecting path, since two would make a cycle. Every off-path bosonic subtree has only external Higgs leaves. A root $W$ or $\chi$ on such a subtree would give odd parity, so its root must be $h$. The far chain couples only through $W$, so it has exactly one vertex; it cannot emit a separate pure-Higgs subtree. Every bosonic vertex on the path keeps the two path edges in the same parity class after accounting for its even side branches. Since the far end is $W$, the path contains only $W/\chi$; near $h$ is emitted only off this path. Distinguishable pairs forbid crossed fermion-line pairings.

With two labelled Higgs leaves the cases are:

| Where the two Higgs bosons are attached | Number |
|---|---:|
| Both on the bridge, two possible orders, two near and two inner field assignments | 8 |
| Both on the bridge via $h^*\to hh$, two near field assignments | 2 |
| $WWhh$ contact (no $\chi Whh$) | 1 |
| One from near line and one on bridge: two labelled assignments, two near-chain orders, two bridge fields | 8 |
| Both directly on massive line: (3!) orders of (W,h_p,h_r) | 6 |
| Both from one line-emitted $h^*\to hh$: two line orders | 2 |

All trees obey (n_3+2n_4=4); no symmetry factor is needed for a fixed externally labelled tree. The factorial for identical legs at a vertex is in the differentiation rule, not an extra diagram multiplicity. An unordered labelled-leaf subset recursion over all neutral cubic and quartic vertices reproduces the table. It roots at the far antifermion and records one canonical signature per graph. This enumeration does not preassign topology classes.

## Exact signed cancellation

Use metric (+---), all incoming momenta and Fourier waves $e^{-ip\cdot x}$. Near bare rules are $\Gamma_W=-ia\gamma^\alpha P_L,\ \Gamma_\chi=b\gamma_5,\ \Gamma_h=-ib$, with $b=\mu/v$. The far current $e^\alpha=\bar v_g\gamma^\alpha u_g$ obeys $K\cdot e=0$. Thus far sewing in every graph gives the same pushed vector $-a e/(K^2-m^2)$. Strip that far propagator normalization and the common near factor $-ia$ in the algebra below; the code includes both. Let $q=p_f+p_{\bar f}$, (q+K+p+r=0), $j^\alpha=\bar v_f\gamma^\alpha P_Lu_f$, and $u_q=q\cdot j=\mu\bar v_f\gamma_5u_f$.

In the bridge-only ordered contribution labelled by $x\in\{p,r\}$, take $Q=q+x$, $D=Q^2-m^2$, $C=Q^2-\xi m^2$. Denote the right subamplitudes with vector and Goldstone near roots by $iT^\alpha$ and $G$. The vector rules give
```text
E·T = 4a²m² [E·e − (1−ξ)(E·Q)(e·Q)/C]/D
      − 2a² [(2Q−q)·E](e·Q)/C,
mG  = 2a²m² [(2q−Q)·e − (1−ξ)((2q−Q)·Q)(e·Q)/C]/D
      + 2a² M²(e·Q)/C.
```
Using $2q\cdot Q-q^2+M^2=Q^2$, each ordering has $q\cdot T-mG=-2a^2 e\cdot Q$. Those terms sum to $-2a^2 e\cdot q$, as (Q_p+Q_r=q-K). The $WWhh$ contact contributes the opposite sign, and the bridge $h^*\to hh$ separately has $q\cdot T=mG$. Thus the full bridge-only subtree has $q\cdot T=mG$ exactly. The delicate point is that the *ordered* pieces do **not** separately obey that Ward relation. With the physical near $W/\chi$ source, the root propagator pair reduces to $[j\cdot T-u_q(q\cdot T)/m^2]/(q^2-m^2)$.

The remaining gauge variation comes from ordered inner exchanges. From the first displayed right-subtree rule, at fixed momenta
```text
E·(Tξ−T1) = 2a² m² (ξ−1)(E·q)(e·Q)/(D C).
```
In the root-sewn sum, the shift from each ordered exchange is therefore $-2a^2(\xi-1)u_q(e\cdot Q)/(D C)$. This is a statement about the variation after simplifying the *summed* root pair, not separate root-ST identities per ordering. The two near-line diagrams for the same (x) have $Q\cdot J_h=\mu P_h+\bar v_f\gamma_5u_f$; the endpoint term follows by commuting the chiral projectors through the two inverse fermion propagators on shell. Their raw $W/\chi$ pair with the far single-Higgs vertex gives
```text
A_mixed = b [J_h·T_far − (Q·J_h)(Q·T_far)/m²]/D
          + b (vbar_f γ5 u_f)(Q·T_far)/(m² C),
Q·T_far = 2 a m e·Q.
```
Its shift is $+2a^2(\xi-1)u_q(e\cdot Q)/(D C)$. It cancels per labelled assignment. The remaining six direct and two split near diagrams are individually $\xi$ independent since only the common far $W$ connects and $K\cdot e=0$. At generic nonsingular kinematics the gauge-independent answer equals its $\xi\to\infty$ limit: the Goldstone propagators decouple and the vector propagators become unitary. No small-width or large-energy limit enters this argument.

## Reproduction and limits

The direct program differentiates the signed neutral Fourier action and contracts all labelled tensor graphs, including fermion propagator orientation. The independent unitary implementation never calls the graph generator or action vertex evaluator; it enumerates four contracted contributions by closed matrix formulas. They share kinematics, on-shell wavefunctions, and a gamma-matrix convention, not gauge assembly. The pre-existing [double-Higgs fixture source](../su2-sixpoint-double-higgs-radiation/check_double_higgs.py) is hash-locked. Five finite common gauges are tested for each seed; fixtures require all relevant absolute propagator denominators greater than 0.22 in recorded units. Exact SymPy checks cover the ordered ST remnant, pure-root-index variation, root sewing and paired mixed variation. The test intentionally removes both the entire mixed class and just the mixed Goldstone graphs: both omissions must introduce a nonzero gauge residual. Native signed [real,imaginary] graph values, fixture momenta and spinors are saved for replay, without phase fits or squaring amplitudes.

```sh
cd research-baselines/su2-sixpoint-connected-tree
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 python check_connected_tree.py --seed 142023 --points 8 --output results-primary.json > summary-primary.txt
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 python check_connected_tree.py --seed 642271 --points 12 --output results-holdout.json > summary-holdout.txt
sha256sum -c SHA256SUMS
```

Source: this branch at [3b7146ab](https://github.com/lm846802297-bit/off-shell-vector-polarization/tree/3b7146abfec3a729aca9431c41e2efc05b591ea9); the rooted-source derivation was separately archived at [5ec05218](https://github.com/lm846802297-bit/off-shell-vector-polarization/tree/5ec052181b83d9b43f1ef4999a1b882538376f6b/research-baselines/su2-fivepoint-yukawa-source). Ke–Liu–Mahbub [arXiv:2611.13096v1 §4.2](https://arxiv.org/html/2611.13096v1#S4.SS2) discuss the distinction between an on-shell construction and off-shell vector–Goldstone constraints. Our check is a Lagrangian tree calculation in one model and channel, not a derivation from their on-shell recursion framework.

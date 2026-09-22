
#import "@preview/touying:0.7.4": *
#import "@preview/shuimu-touying-zen:0.1.0": *


#show: group-meeting-theme.with(
  cover-logo-name: "phys-logo.svg",
  config-common(handout: sys.inputs.at("handout", default: "false") == "true"),
  config-info(
    author: "Chenyang Wang",
    institution: "Tsinghua University \n Email: cy-wang21@mails.tsinghua.edu.cn",
    subtitle: "",
    date: datetime.today(),
    title: "Numerical computation of 2D \n generalized Brillouin zones"
  )
)

#import "@preview/lovelace:0.3.1": *

#let my-pseudo = pseudocode-list.with(line-gap:0.8em, booktabs: true)

#let mathbf(x) = $bold(upright(#x))$
#let rme = $upright(e)$
#let rmi = $upright(i)$
#let reference-item(
  author, journal, volume-pages, year
) = text(18pt)[
  #author, #emph(journal), #volume-pages (#year)
]

// Keep the existing theme and reveal order. Use --input handout=true for static slides.


#set list(spacing: 1em)
#set enum(spacing: 1em)
#set par(spacing: 1.2em)
#show math.equation.where(block: true): set block(above: 0.5em, below: 0.5em)

#title-slide()

#outline-slide()

= Introduction to 2D generalized Brillouin zones

== Non-Bloch band theory
#grid(columns: (1fr,1fr), gutter: 32pt,
[
  *Bloch band theory*
  #pause
  - Hermitian periodic bulk
    \ $psi(x+a)=exp(rmi k a)psi(x)$
  #pause
  - Bulk spectrum in the OBC limit
    \ Bloch energy bands $E_n(k)$
  #pause
  - Bulk OBC states are superpositions of Bloch waves $psi_(k,n)(x)$.
  #align(center)[#image("Figures/schematic-Bloch.pdf", height: 135pt)]
],
[
  #meanwhile
  *Non-Bloch band theory*
  #pause
  #pause
  #pause
  #pause
  - #text(red)[Non-Hermitian] periodic bulk
    \ $psi(x+a)=beta psi(x), beta in CC, beta eq.not 0$
  #pause
  - Bulk spectrum in the OBC limit
    \ Non-Bloch energy bands $E_n(beta)$
  - Bulk OBC states are superpositions of non-Bloch waves $psi_(beta,n)(x)$.
  #align(center)[#image("Figures/schematic-non-Bloch.pdf", height: 135pt)]
])

== Generalized Brillouin zone

- For an $n$-dimensional lattice, extending $mathbf(k) in RR^(n) -> bold(beta) in CC^(n)$ *doubles the dimensionality* of momentum space.
  #pause
- $n$ real-valued constraints on $bold(beta)$ are required to recover the dimensionality.
  #pause
- 1D generalized Brillouin zone (GBZ): equal decay rates of the boundary-selected waves.
  #pause
  - Write $det (E - h (beta)) = a_(- M) beta^(- M) + a_(- M + 1) beta^(- M + 1) + dots.c + a_(N) beta^(N)$
  #pause
  - $det ( E - h ( beta ) ) = 0 =>$ sort the roots by modulus: $|beta^((1))| <=|beta^((2))| <= dots.c <= |beta^((M+N))|$
  #pause
  - GBZ constraint: $lr(|beta^((M))|)=lr(|beta^((M+1))|)$

  #columns(2)[
    #image("Figures/GBZ-YaoWang.png", height: 100pt)
    #colbreak()
    #align(horizon)[
      - #reference-item(
        "Yao, S. & Wang, Z.",
        "\n Phys. Rev. Lett.",
        "121, 086803", 2018
      )
      - #reference-item(
        "Yokomizo, K. & Murakami, S.", "\nPhys. Rev. Lett.", "123, 066404", "2019"
      )
    ]
  ]

== Amoeba formulation
- In two and higher dimensions, the zeros of $det [E - h (bold(beta))]$ are not discrete.
  \ A global modulus ordering no longer yields a direct analogue of the 1D condition.
#pause
- Amoeba formulation: $ E in sigma_("Amoeba") <=> "Absence of a \"central hole\" in" cal(A)(E) $

#grid(
  columns: (2fr, 1fr),
  [
    #align(horizon)[
      - $cal(A) (E) := {bold(mu) in RR^(n) | exists bold(theta) , det [E - h (rme^(bold(mu) + rmi bold(theta)))] = 0}$
      #pause
      - Ronkin function: $ R_E (bold(mu))=integral_(TT^n)(frac(dif bold(theta), 2 pi))^(n)ln lr(|det[E-h("e"^(bold(mu) + "i"bold(theta)))]|) $
      - Central hole: an open complement component where $nabla R_E = 0$.
    ]
  ],
  [
    #meanwhile
    #pause
    #align(center)[
      #image("Figures/WangAmoeba.png")

      #reference-item("Wang, H.-Y., Song, F. & Wang, Z.", "Phys. Rev. X","14, 021011", "2024")
    ]
  ]
)

== Strip-GBZ formulation
#columns(2)[
*Geometry-dependent skin effect*
- In $n >= 2$ dimensions, the OBC spectrum and decay factors depend on boundary geometry.
#align(center)[#image("Figures/GDSE-20250402.pdf", width: 100%)
]

#colbreak()
*Strip-GBZ formulation*
- OBC spectrum in the strip limit ($L_1 >> L_2$).
#align(center)[
  #image("Figures/SGBZ-formulation-20260325.pdf", width:80%)
]
]

#columns(2)[
 #reference-item("Zhang, K., et al.", "Nat. Commun.", "13, 2496", "2022")
 #colbreak()
 #align(center)[
  #reference-item("Xiong, Y., et al.", "arXiv", "2407.01296", "2024")
  \ #reference-item("Wang, C., et al.", "arXiv", "2506.22743", "2025")]
]

#slide()[
  #grid(
    columns: (1.5fr, 1fr),
    align: (top+left, center+horizon),
    [
      #pause
      - Minor-axis log-radius function: $rho_(2 , 0) (theta_(1) ; E , mu_(1))$

        Periodic, with $mu_2^((M_2))<= rho_(2,0) <= mu_2^((M_2+1))$
      #pause
      - Base manifold:  $ X (E , mu_(1)) := & lr({(beta_(1) , beta_(2)) | theta_(1) , theta_(2) in [- pi , pi] \
 & beta_(1) = "e"^(mu_(1) + "i" theta_(1)) , beta_(2) = "e"^(rho_(2 , 0) (theta_(1)) + "i" theta_(2))}) $
      #pause
      - Major-axis winding number $w_1(theta_2;E,mu_1)$:
        \ Winding of $det(E-h)$ along $theta_1$ at fixed $theta_2$.
      #pause
      - Average winding number:
        $ W(E,mu_1) := frac(1, 2 pi) integral_(-pi)^(pi)w_1(theta_2;E,mu_1) dif theta_2 $
    ],
    [
      #meanwhile
        #image("Figures/SGBZ-formulation-20250818.pdf", width:100%)
    ]
  )
]

== Relation between the two frameworks
- Derivative of $R_E (bold(mu))$ $<=>$ average winding number on the torus at fixed $bold(mu)$
  $ frac(partial R_E (bold(mu)), partial mu_j) = integral_(TT^(n-1)) frac(dif theta_1 dots.c hat(dif theta_j) dots.c dif theta_n, (2 pi)^(n-1)) u_j (theta_1, ..., hat(theta_j), ..., theta_n) $
  where
  #v(-1em)
  $ u_(j) (theta_(1) , ... , hat(theta_(j)) , ... , theta_(n)) := integral_(0)^(2 pi) frac(dif theta_(j), 2 pi rmi) frac(partial_(theta_(j)) det (E - h (rme^(bold(mu) + rmi bold(theta)))), det (E - h (rme^(bold(mu) + rmi bold(theta))))) $

#align(center)[
  #image("Figures/Amoeba-20260207.pdf", width: 90%)
]


= Computing average winding via polynomial root-solving

== Winding number as a homomorphism $H_1 -> ZZ$
- $det(E-h(bold(beta)))$: a *continuous* map from the base manifold $bold(beta) in X$ to the complex plane.
- Projection onto the unit circle: $P:CC \\{0} -> S^1, z mapsto z \/|z|$

#grid(
  columns: (1fr, 1fr),
[
- Composition:
  $ X \\ {"Zeros"} stretch(->)^(det [E - h (bold(beta))]) CC \\ {0} stretch(->)^(P) S^(1) $
- Winding number: homomorphism
  $ w_("loop"): pi_1(X\\{"Zeros"}) --> pi_1(S^1)=ZZ $
- Since $ZZ$ is Abelian,
  $ w_("loop"): H_1(X\\{"Zeros"}) --> ZZ $

],
[
  #align(center + horizon)[#image("Figures/winding-general-loop-20250527.pdf", width: 12cm)]
  - Strip GBZ: $X=X(E,mu_1)$
  - Amoebic GBZ: $X = X (bold(mu)) := {(rme^(mu_(1) + rmi theta_(1)) , rme^(mu_(2) + rmi theta_(2))) | theta_(1) , theta_(2) in [- pi , pi]}$
])

== Computation of $nabla R_E$ (Amoeba)
- Argument principle: $ integral.cont_(cal(C)) frac(dif z, 2 pi rmi) frac(f'(z), f(z)) = N_("zeros") - N_("poles") $
  Zeros and poles inside $cal(C)$ count with multiplicity. The contour avoids both.
- Let $-M_j$ and $N_j$ be the lowest and highest powers of $beta_j$ in $det(E-h)$, for $j=1,2$.
- Loop winding number:
  - $u_(1) (theta_(2) ; mu_(1) , mu_(2)) = ("# of" beta_(1) "roots inside" rme^(mu_(1))) - M_(1)$
  - $u_(2) (theta_(1) ; mu_(1) , mu_(2)) = ("# of" beta_(2) "roots inside" rme^(mu_(2))) - M_(2)$

#slide()[
#grid(columns: (2fr, 1fr),
[
	#my-pseudo(title:[*Algorithm*: Compute $partial_(2) R_(E) (mu_(1) , mu_(2))$])[
		+ Compute tracks of roots $(E , mu_(1) , theta_(1)) mapsto beta_(2)^(( j ))$
		+ Collect the intersection points of $ln |beta_(2)^((j))| = mu_(2)$
		+ *if* intersection points are discrete:
			+ Partition the $theta_(1)$ circle at the crossing angles
			+ *For each* interval:
				+ Choose the interval midpoint $theta_(1)$
				+ Count roots with $ln|beta_2|<mu_2$: $u_2=N-M_2$
			+ $partial_(2) R_(E) = sum_("intervals") u_(2) times ( "interval width" ) \/ 2 pi$
		+ *else*: $partial_(2) R_(E)$ is ill-defined.
	]

],
[#image("Figures/Ronkin-calc-20260914.pdf", width: 100%)
  - $u_(1) (theta_(2) ; mu_(1) , mu_(2)) =$ \ $("# of " ln | beta_(1) | < mu_(1)) - M_(1)$
  - $u_(2) (theta_(1) ; mu_(1) , mu_(2)) =$ \ $("# of " ln | beta_(2) | < mu_(2)) - M_(2)$
	#pause
	- #text(blue)[ $partial_(1) R_(E)$ can be computed with the same procedure.]
]

)
]

== Discrete vs. continuum
#columns(2)[
  *Discrete zeros (point subsets)*
	#image("Figures/discrete-subsets.pdf", width: 80%)
	- Singular loop parameters form a set of measure zero.
	- $integral_(- pi)^(pi) u_(2) (theta_(1)) dif theta_(1)$ is #text(blue)[well defined].

	#colbreak()

	*Continuum zeros (line subsets)*
	#image("Figures/continuum-subsets.pdf", width: 80%)
	- Singular loop parameters occupy an interval.
	- $integral_(- pi)^(pi) u_(2) (theta_(1)) dif theta_(1)$ is #text(red)[ill-defined].
]

#slide()[
  #columns(2)[
  *Discrete zeros (point subsets)*
	#image("Figures/discrete-partial-R.pdf", width: 80%)

	#colbreak()

	*Continuum zeros (line subsets)*
	#image("Figures/continuum-partial-R.pdf", width: 80%)
]]


== Nested bisection for the Ronkin minimum
- At fixed $E$, GBZ subsets satisfy $f(E,beta_1,beta_2)=0$ and $ln |bold(beta)| = bold(mu)_min$.
- Partial minimization preserves convexity: $F_E (mu_1) := min_(mu_2) R_E (mu_1,mu_2)$ is convex.

#grid(columns: (1fr, 1fr), gutter: 28pt,
[
  *Outer search in $mu_1$*
  #my-pseudo(line-gap: 0.8em)[
    + Bracket the sign change in $partial_1 R_E$.
    + Choose the midpoint $mu_1$.
    + Track $beta_2^((j))(theta_1)$ once.
    + Run the inner search for $mu_2$.
    + Evaluate $partial_1 R_E$ at that result.
    + Update the $mu_1$ bracket.
    + Stop at zero winding or a resolved continuum minimum.
  ]
],
[
  *Inner search in $mu_2$*
  #my-pseudo(line-gap: 0.8em)[
    + Reuse the tracks at fixed $(E,mu_1)$.
    + Test the gap and flat-track levels.
    + Bracket the sign change in $partial_2 R_E$.
    + Locate crossings of $ln|beta_2|=mu_2$.
    + Evaluate $partial_2 R_E$ by root counts.
    + Bisect until $|partial_2 R_E| < epsilon_w$.
    + Retain refined roots for later levels.
  ]
])

== Continuum zeros and one-sided winding limits
- A flat root track at $ln|beta_2|=mu_(2,c)$ makes the usual loop winding singular over an interval of $theta_1$.
- Probe the two sides of the level:
  $ w_2^- = partial_2 R_E (mu_1,mu_(2,c)-epsilon), quad
    w_2^+ = partial_2 R_E (mu_1,mu_(2,c)+epsilon). $
- If $w_2^- < 0 < w_2^+$, the level minimizes the convex function in $mu_2$.
  The minimum can occur at a kink where an ordinary derivative does not exist.
- Resolve the outer search with analogous $mu_1$ probes at fixed $mu_(2,c)$.
  Opposite signs identify a continuum minimum and produce line subsets.
- If the signs do not straddle zero, continue on the appropriate side.
  A zero-winding plateau with no zeros represents an exterior energy.

= Pseudo-arclength continuation of zeros

== Root tracks at fixed energy and decay rate
- Fix $(E,mu_1)$ and vary $theta_1$, with $beta_1=exp(mu_1+rmi theta_1)$.
  Solve $f(E,beta_1,beta_2)=0$ for all $K=M_2+N_2$ roots.
- A uniform $theta_1$ mesh can underresolve rapid root motion near a multiple root.
  Sorting roots independently at each step also loses their track identities.
- The continuation state is
  $ bold(y)=(theta_1, ln beta_2^((1)), dots.c, ln beta_2^((K))). $
- Adapt the angular step to motion in this state space. At each new angle,
  solve the polynomial again and match the roots to their predicted positions.

== Tangent and arclength step
#show math.equation.where(block: true): set block(above: 0.5em, below: 0.5em)
Implicit differentiation along a regular root gives
$ partial_(beta_1)f dot rmi beta_1 + partial_(beta_2)f dot frac(dif beta_2,dif theta_1)=0. $
Hence the logarithmic tangent of track $j$ is
$ V_j := frac(dif ln beta_2^((j)),dif theta_1)
  = -rmi frac(beta_1, beta_2^((j)))
    frac(partial_(beta_1)f, partial_(beta_2)f), quad
  bold(V)=(1,V_1,dots.c,V_K). $
The proposed angular increment is
$ Delta theta_1=frac(h,norm(bold(V))_2), quad
  norm(bold(V))_2=sqrt(1+sum_j |V_j|^2). $
- Fast root motion increases $norm(bold(V))_2$ and reduces $Delta theta_1$.
- $h$ controls the approximate arclength step and adapts to prediction error.
- Undefined tangents of zero/infinite padding roots are excluded from the norm.

== Prediction, polynomial solve, and track matching
#my-pseudo(line-gap: 0.6em)[
  + *Predict*: $hat(beta)_2^((j))=beta_2^((j)) exp(V_j Delta theta_1)$.
  + *Solve*: compute all roots $tilde(beta)_2^((k))$ at the trial angle.
  + *Match*: find a one-to-one assignment $pi$ minimizing
    $ sum_j d_("chord")(hat(beta)_2^((j)),tilde(beta)_2^((pi(j)))). $
  + *Accept or retry*: compare matched roots with the predictions.
]

#v(-.4em)
- Hungarian matching uses distance on the Riemann sphere:
  $ d_("chord")(z,w)=frac(2|z-w|,sqrt(1+|z|^2)sqrt(1+|w|^2)). $
- The sphere metric handles roots near zero and infinity.
  A prediction anchor reduces track swaps at close approaches.
- Each regular root sample comes from a polynomial solve.

== Adaptive step-size control
Use the largest matched prediction error:
$ e=frac(max_j d_("chord")(hat(beta)_2^((j)),tilde(beta)_2^((pi(j)))),
  "atol"+"rtol" dot limits("median")_(k:"finite") |tilde(beta)_2^((k))|). $

#v(.4em)

#grid(columns: (1fr,1fr), gutter: 30pt,
[
  *Accepted step: $e<1$*
  - Store the solved roots in track order.
  - Propose $h_("next") = h min(10,0.9 e^(-1/2))$.
  - For $e=0$, use the maximum growth factor.
],
[
  *Rejected step: $e>=1$*
  - Reduce $h$ and solve again.
  - Use $h <- h max(0.2,0.9 e^(-1/2))$.
  - Stop after the retry budget or minimum step is reached.
])
#v(0.5em)
Defaults: $"atol"=10^(-12)$, $"rtol"=10^(-3)$, $h_("max")=0.5$.
An accepted step cannot grow immediately after a rejection.

== Multiple roots split the continuation tracks
#table(
  columns: (1fr, 2fr), inset: 9pt, stroke: 0.5pt + luma(75%),
  table.header([*Trigger*], [*Candidate condition*]),
  [Step collapse], [$Delta theta_1 < 10^(-10)$ near a singular implicit derivative.],
  [Pairwise approach], [For a nearby root pair, $dif |beta_i-beta_j|^2 \/ dif theta_1$ changes from negative to positive.],
)
- Refine a point candidate by solving $f=0$ and $partial_(beta_2)f=0$.
  For an interval candidate, first locate its distance-derivative zero.
- A close approach is only a candidate. Confirm a finite-root cluster before splitting the segment.
- Store the multiple root as a shared endpoint, restart beyond it, and track the next segment.
- At $theta_1=2pi$, retain the permutation relating the closing roots to $theta_1=0$.

== Root tracking shared by both GBZ solvers
- `ZeroManager` returns tracked roots, analytic tangents, and multiple-root boundaries on an adaptive mesh.
- *SGBZ*: refine equal-modulus events, build $mu_(2,"mid")(theta_1)$, and extract boundary points or line intervals.
- *Amoeba*: locate crossings of $ln|beta_2|=mu_2$ and reuse the same tracks throughout the inner $mu_2$ search.
- Cubic Hermite interpolation supplies matching anchors and the SGBZ loop path.
  Singular endpoint derivatives require linear fallback.
- Step control resolves root motion. Additional crossing and extremum refinement resolves the level sets used by each solver.

= Benchmark: 2D complex Hatano-Nelson model

== Hopping parameters and geometry
$ h(beta_x,beta_y)=J_(x 1) beta_x^(-1)+J_(x 2) beta_x
                 +J_(y 1) beta_y^(-1)+J_(y 2) beta_y. $
#grid(columns: (1fr,1.2fr), gutter: 28pt,
[
  For $alpha=x,y$, write the nonzero hoppings as
  $ J_(alpha 1)=exp(gamma_alpha+rmi delta_alpha)J_alpha, $
  $ J_(alpha 2)=exp(-gamma_alpha+rmi delta_alpha)J_alpha^*. $
  - $gamma_alpha$: nonreciprocal decay rate.
  - $delta_alpha$: common hopping phase.
  - Compare amoeba with $x$-, $y$-, and $[11]$-strip GBZs.
],
[
  #align(center + horizon)[#image("../../application/Figures/HN-example-20250524.pdf", width: 100%)]
])
#v(0.5em)
#reference-item("Wang, C., et al.", "arXiv", "2506.22743v3, Eqs. (S3.16), (S3.30)", "2025")

== Closed-form GBZ constraints
#show math.equation.where(block: true): set block(above: 0.5em, below: 0.5em)
*Cartesian coordinates: amoeba, $x$ strip, and $y$ strip*
$ |beta_x|=exp(gamma_x), quad |beta_y|=exp(gamma_y). $
*Diagonal strip: $(beta_1,beta_2)=(beta_x beta_y,beta_y)$*
$ |beta_1|=exp(gamma_x+gamma_y), quad
  A=frac(J_(x 1),beta_1)+J_(y 2), quad B=J_(x 2)beta_1+J_(y 1). $
At fixed $beta_1$, the characteristic equation becomes
$ A beta_2^2-E beta_2+B=0, quad |beta_2|=sqrt(lr(|B/A|)). $
- The product of the two roots is $B/A$. Equal moduli give the transverse GBZ radius.
- The $[11]$ radius generally varies with $theta_1$. These formulas apply where $A$ and $B$ are nonzero.

== Accuracy measures and test cases
For every returned point and every stored line sample, compare
$ epsilon_mu=max_(j=1,2) max_n |ln|beta_(j,n)|-mu_j^"exact"(theta_(1,n))|, $
$ epsilon_E=max_n frac(|E-h(beta_(x,n),beta_(y,n))|,
 |E|+|J_(x 1)/beta_(x,n)|+|J_(x 2)beta_(x,n)|+|J_(y 1)/beta_(y,n)|+|J_(y 2)beta_(y,n)|). $
#table(
  columns: (auto,1fr,auto), inset: 8pt, stroke: 0.5pt + luma(75%),
  table.header([*Case*], [*$(J_(x 1),J_(x 2),J_(y 1),J_(y 2))$*], [*$E$*]),
  [Point subsets], [$(1+rmi,1.5+1.2rmi,-1+rmi,-1.2-0.5rmi)$], [$1+rmi$],
  [Line subsets], [$(1,1.5,-1,-1.2)$], [$1$],
)
- Acceptance thresholds: $epsilon_mu <= 10^(-5)$ and $epsilon_E <= 10^(-7)$.
- Check spectral membership independently, including successful empty results.

== Numerical errors for all four GBZs
#set par(spacing: 0.6em)
#set text(size: 20pt)
// Current-checkout results from application/benchmark-2D-Hatano-Nelson.py.
// Recomputed 2026-09-22 with NumPy 1.26.4 and SciPy 1.11.4.
#table(
  columns: (auto,auto,auto,auto,auto), inset: (x: 13pt, y: 6pt),
  stroke: 0.5pt + luma(75%), align: (left,left,right,right,right),
  table.header([*Case*], [*GBZ*], [*Samples*], [*$epsilon_mu$*], [*$epsilon_E$*]),
  [Points], [Amoeba], [4], [$9.52 times 10^(-9)$], [$2.01 times 10^(-16)$],
  [Points], [$x$ strip], [4], [$2.00 times 10^(-9)$], [$1.00 times 10^(-16)$],
  [Points], [$y$ strip], [4], [$1.44 times 10^(-9)$], [$1.00 times 10^(-16)$],
  [Points], [$[11]$ strip], [4], [$1.44 times 10^(-9)$], [$2.10 times 10^(-16)$],
  [Lines], [Amoeba], [354], [$5.97 times 10^(-7)$], [$1.32 times 10^(-10)$],
  [Lines], [$x$ strip], [350], [$4.22 times 10^(-7)$], [$1.32 times 10^(-10)$],
  [Lines], [$y$ strip], [346], [$5.80 times 10^(-7)$], [$1.28 times 10^(-10)$],
  [Lines], [$[11]$ strip], [380], [$3.98 times 10^(-7)$], [$5.51 times 10^(-11)$],
)
#v(0.5em)
All eight cases pass. Each point case contains four points and each line case contains two line subsets.

#text(size: 17pt)[NumPy backend, default solver settings. Sample accuracy does not by itself prove that every branch was found.]

== Spectral membership beyond nonempty examples
#set par(spacing: 0.7em)
#set list(spacing: 0.8em)
Let $u=2exp(rmi delta_x)|J_x|$ and $v=2exp(rmi delta_y)|J_y|$.
- Cartesian spectra: $E in {u s+v t: s,t in [-1,1]}$.
- Diagonal strip: $E^2 in "conv"{0,(u+v)^2,(u-v)^2}$.
#table(
  columns: (auto,auto,1fr,1fr), inset: 9pt, stroke: 0.5pt + luma(75%),
  table.header([*Hoppings*], [*$E$*], [*Amoeba / $x$ / $y$*], [*$[11]$ strip*]),
  [Complex case], [$6+6rmi$], [Outside], [Outside],
  [Real case], [$6$], [Outside], [Outside],
  [Real case], [$rmi$], [Outside], [Outside],
  [Complex case], [$2+2rmi$], [Inside], [Outside],
)
Sixteen solver comparisons agree with the independent analytic classification:
13 empty results and 3 nonempty results. Failed solves do not count as exterior energies.

== Energy-grid benchmark
#set par(spacing: 0.55em)
// Plotted from the four application/data/HN2D-20260921T144330717373Z-*.pkl files.
#align(center)[#image("Figures/hn-benchmark-spectra.pdf", height: 250pt)]
- All four methods evaluate the same $100 times 100$ grid, with zero failed results.
- Amoeba, $x$ strip, and $y$ strip share all 4,832 in-spectrum grid points.
  The $[11]$ strip contains 1,108 points on the same grid.
#text(size: 17pt)[Complex hopping case. Dots are computed grid points. Finite sampling does not certify arbitrarily narrow spectral features.]

// Application slides are reserved for the next part of the talk.
// = Applications
// == Geometry-dependent skin effect
// == Chern number calculations

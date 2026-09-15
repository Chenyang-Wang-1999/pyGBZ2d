
#import "@preview/touying:0.7.4": *
#import "@preview/shuimu-Touying-zen:0.1.0": *


#show: group-meeting-theme.with(
  cover-logo-name: "phys-logo.svg",
  config-info(
    author: "Chenyang Wang",
    institution: "Affil: Tsinghua University \n Email: cy-wang21@mails.tsinghua.edu.cn",
    subtitle: "",
    date: datetime.today(),
    title: "Numerical calculation of 2D \n generalized Brillouin zone"
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

// Rendering speed-up


#set list(spacing: 1.2em)
#set enum(spacing: 1.2em)
#set par(spacing: 1.5em)
#show math.equation.where(block: true): set block(above: 1em, below: 1em)

#title-slide()

#outline-slide()

= Introduction to 2D generalized Brillouin zone

== Non-Bloch band theory

#columns(2)[
  *Bloch band theory*
  #pause
  - Hermiticity + Periodicity 
    \ $=> psi (x + a) = rme^(rmi k a) psi (x)$
  #pause
  - OBC spectrum  
    \ $->$ Bloch energy bands $E_n (k)$
  #pause
  - OBC eigenstates 
    \ $->$ Bloch wavefunctions $psi_(k,n) (x)$

  #move(dx: 20pt)[#image("Figures/schematic-Bloch.pdf", width: 75%)]

  #colbreak()
  #meanwhile
  *Non-Bloch band theory*
  #pause
  #pause
  #pause
  #pause
  - #text(red)[Non-Hermiticity] +  Periodicity
    \ $=> psi (x + a) = beta^(a) psi (x) 、 , beta in CC$

  #pause
  - OBC spectrum
    \ $->$ Non-Bloch energy bands $E_n (beta)$
  - OBC eigenstates
    \ $->$ Non-Bloch wavefunctions $psi_(beta , n) (x)$
  
  #align(center)[#image("Figures/schematic-non-Bloch.pdf", width:75%)]

]

== Generalized Brillouin zone

- For $n$D lattice, extension $mathbf(k) in RR^(n) -> bold(beta) in CC^(n)$ *doubles the dimensionality* of momentum space.
  #pause
- $n$ real-valued constraints on $bold(beta)$ are required to recover the dimensionality.
  #pause
- 1D generalized Brillouin zone (GBZ): "coherence" of non-Bloch waves.
  #pause
  - Assuming $det (E - h (beta)) = a_(- M) beta^(- M) + a_(- M + 1) beta^(- M + 1) + dots.c + a_(N) beta^(N)$
  #pause
  - $det ( E - h ( beta ) ) = 0 =>$ roots sorted in order of $|beta^((1))| <=|beta^((2))| <= dots.c <= |beta^((M+N))|$
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
- For 2 and higher dimensions, zeros of $det [E - h (bold(beta))]$ are not discrete.
  \ The constraint $| beta^((M)) | = | beta^((M + 1)) |$ lacks higher-dimensional generalization.
#pause
- Amoeba formulation:$ E in sigma_("GBZ") <=> "Absence of a \"central hole\" in" cal(A)(E) $

#grid(
  columns: (2fr, 1fr),
  [
    #align(horizon)[
      - $cal(A) (E) := {bold(mu) in RR^(n) | exists bold(theta) , det [E - h (rme^(mu + "i" bold(theta)))] = 0}$
      #pause
      - Ronkin function: $ R_E (bold(mu))=integral_(TT^n)(frac(dif bold(theta), 2 pi))^(n)ln lr(|det[E-h("e"^(bold(mu) + "i"bold(theta)))]|) $
      - Central hole: constant plateau of $R_E (bold(mu))$
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
- $>=2$D: OBC spectrum & $bold(beta)$ depends on the geometry of boundaries.
#align(center)[#image("Figures/GDSE-20250402.pdf", width: 100%)
]

#colbreak()
*Strip-GBZ formulation*
- Limit case of OBC spectrum in "strip geometries" ($L_1 >> L_2$)
#align(center)[
  #image("Figures/SGBZ-formulation-20260325.pdf", width:80%)
]
]

#columns(2)[
 #reference-item("Zhang, K., et al.", "Nat. Commun.", "13, 2496", "2022")
 #colbreak()
 #align(center)[
  #reference-item("Xiong, Y., et al.", "arXiv", "2407.01296", "2024") 
  \ #reference-item("Wang C. et al.", "arXiv", "2506.22743", "2025")]
]

#slide()[
  #grid(
    columns: (1.5fr, 1fr),
    align: (top+left, center+horizon),
    [
      #pause
      - Minor-axis radius function: $rho_(2 , 0) (theta_(1) ; E , mu_(1))$

        Periodic function s.t. $mu_2^((M_2))<= rho_(2,0) <= mu_2^((M_2+1))$
      #pause
      - Base manifold:  $ X (E , mu_(1)) := & lr({(beta_(1) , beta_(2)) | theta_(1) , theta_(2) in [- pi , pi] \
 & beta_(1) = "e"^(mu_(1) + "i" theta_(1)) , beta_(2) = "e"^(rho_(2 , 0) (theta_(1)) + "i" theta_(2))}) $
      #pause
      - Major-axis winding number $w_1(theta_2;E,mu_1)$:
        \ Winding number on $X(E,mu_1)$ around $theta_1$ axis
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
- Derivative of $R_E (bold(mu))$ $<=>$ average winding number on the plane with constant $bold(mu)$
  $ frac(partial R_E (bold(mu)), partial mu_j) = integral_(TT^(n-1)) frac(dif theta_1 dots.c hat(dif theta_j) dots.c dif theta_n, (2 pi)^(n-1)) u_j (theta_1, ..., hat(theta_j), ..., theta_n) $
  where
  #v(-1em)
  $ u_(j) (theta_(1) , ... , hat(theta_(j)) , ... , theta_(n)) := integral_(0)^(2 pi) frac(dif theta_(j), 2 pi rmi) frac(partial_(theta_(j)) det (E - h (rme^(bold(mu) + rmi bold(theta)))), det (E - h (rme^(bold(mu) + rmi bold(theta))))) $
  #v(-1em)
#align(center)[
  #image("Figures/Amoeba-20260207.pdf")
]


= Computing average winding via polynomial root-solving

== Winding number as homomorphism $H_1 -> ZZ$
- $det(E-h(bold(beta)))$: a *continuous* map from the base manifold $bold(beta) in X$ to a complex plane.
- Projection to unit circle: $P:CC \\{0} -> S^1, z mapsto z \/|z|$

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
- Cauchy argument principle: $ integral.cont_(cal(C)) frac(dif z, 2 pi) frac(partial f(z), f(z)) = N_("zeros") - N_("poles") $
  $N_"zeros"$: \# of zeros enclosed by $cal(C)$,
  $N_"poles"$: total degree of poles enclosed by $cal(C)$
- Assuming $beta_j^(-M_j)$ and $beta_j^(N_j)$ are the terms with lowest and highest degrees \ of $beta_j, j=1,2$ in $det(E-h(beta_1,beta_2))$ 
- Loop winding number:
  - $u_(1) (theta_(2) ; mu_(1) , mu_(2)) = ("# of" beta_(1) "zeros below" rme^(mu_(1))) - M_(1)$
  - $u_(2) (theta_(1) ; mu_(1) , mu_(2)) = ("# of" beta_(2) "zeros below" rme^(mu_(2))) - M_(2)$

#slide()[
#grid(columns: (2fr, 1fr),
[
	#my-pseudo(title:[*Algorithm*: Compute $partial_(2) R_(E) (mu_(1) , mu_(2))$])[
		+ Compute tracks of roots $(E , mu_(1) , theta_(1)) mapsto beta_(2)^(( j ))$
		+ Collect the intersection points of $ln ( beta_(2)^((j)) ) = mu_(2)$
		+ *if* intersection points are discrete:
			+ Split $theta_(2)$ axis into intervals by zeros
			+ *For each* interval: 
				+ Find a representative $theta_(1)$
				+ Count zeros below $mu_(2)$, $u_(2) (theta_(1)) = N_("zeros") - M_(2)$
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

== Bisect to $mu_(min)$
- Subset of amoebic GBZ with eigenenergy $E$ $<=>$ zeros satisfying $ln (bold(beta)) = bold(mu)_(min)$ .



= Benchmark: 2D complex Hatano-Nelson model

= Applications

== Chern number calculations

== Dimensional surprise

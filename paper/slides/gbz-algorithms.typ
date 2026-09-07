#import "@preview/touying:0.7.4": *
#import themes.university: title-slide
#import "theme.typ": *
#show: academic-theme
#set document(title: "From physical GBZ formulations to polynomial roots", author: "pyGBZ2d research project", date: none)

#title-slide()
#speaker-note[0:30. The project translates two physical GBZ prescriptions into polynomial roots, crossing geometry, and scalar constraints. This is a talk about numerical design, not a claim of a new GBZ theory.]

= Physical picture

== Complex momentum describes skin localization

#eq[$ psi_(n_1,n_2) prop beta_1^(n_1) beta_2^(n_2), quad beta_j=e^(mu_j+i theta_j) $]
#two([
  #label[Phase and envelope]
  $theta_j$ controls spatial oscillation. $mu_j$ controls exponential growth or decay along a lattice direction.
  #v(16pt)
  A GBZ selects the complex momenta compatible with the boundary problem.
],[
  #label[The 2D question]
  Which pairs of localization exponents and phases survive when boundaries are opened?
  #v(16pt)
  Amoebic and strip constructions organize this selection in different ways.
])
#v(1fr)
#source[Generalized Bloch-wave picture; C. Wang et al., arXiv:2506.22743v2, Sec. II]
#speaker-note[1:00. Start with the wave envelope, since the audience already knows the NHSE. A single generalized Bloch wave is a building block; boundary eigenstates are superpositions. Do not suggest one pair of exponents universally describes every finite geometry.]

== Amoebic GBZ: select localization in modulus space

#note[Write $f(E,beta_1,beta_2)=det[E I-h(beta_1,beta_2)]$.]
#eq[$ cal(A)_E = {(ln abs(beta_1),ln abs(beta_2)): f(E,beta_1,beta_2)=0} $]
#two([
  #label[Project away the phases]
  The amoeba records where a reference energy can be supported by complex Bloch waves with given envelopes.
  #v(14pt)
  The Ronkin minimum selects the logarithmic radii; the GBZ retains the corresponding zeros and their phases.
],[
  #definition-block([Central-hole criterion],[A central zero-free hole gives a flat Ronkin minimum. Its disappearance identifies the amoeba spectral support.])
  #v(12pt)
  #note[Finding a minimum and finding a GBZ point are distinct questions.]
])
#v(1fr)
#source[H.-Y. Wang, F. Song, and Z. Wang, Phys. Rev. X 14, 021011 (2024)]
#speaker-note[1:20. Define f as det(EI-h). The projection loses phase information, so the GBZ needs a lift back to complex Bloch factors. A zero-gradient hole is not a set of supported waves: it is zero-free. The spectral criterion concerns the central hole, not every complement component. At singular minima, use the convex minimum condition rather than presuming a smooth gradient.]

== SGBZ: resolve the boundary problem by a strip

#two([
  #image("Figures/strip-geometry.svg", width:100%)
],[
  #label[Keep a preferred direction]
  Open the transverse direction and treat the strip as a quasi-1D system along its major axis.
  #v(12pt)
  Its thermodynamic winding selects the longitudinal radius. Transverse modes obey the 1D GBZ condition.
  #v(12pt)
  #note[The chosen major direction carries physical information about the boundary geometry.]
], widths:(1.05fr,1fr))
#v(1fr)
#source[C. Wang et al., arXiv:2506.22743v2, Sec. II]
#speaker-note[1:20. Explain the order of the strip construction. The transverse OBC problem provides modes of an effective quasi-1D strip. The strip winding is the width-normalized thermodynamic quantity, not a finite-strip winding without normalization. We will replace this construction by root constraints.]

== Two physical prescriptions, one numerical objective

#two([
  #label[Amoebic GBZ]
  Select radii through a convex spectral potential in the two-dimensional modulus plane.
  #v(18pt)
  #label[SGBZ]
  Select radii through a direction-dependent strip boundary problem.
],[
  #definition-block([Project objective],[Translate each prescription into roots of a one-variable polynomial, while retaining its own physical selection rule.])
  #v(12pt)
  Compare spectra and localization on a common representation.
])
#v(1fr)
#source[The SGBZ theory gives spectral inclusion within the amoeba spectrum; it does not equate the two prescriptions for every geometry.]
#speaker-note[1:00. The common representation does not imply identical GBZs. Distinguish our computational bridge from the theoretical relation between the formulations. The rest of the talk explains how the integral and strip definitions become root problems.]

= The bridge to polynomial roots

== Fix one Bloch factor, solve for the other

#eq[$ f(E,beta_1,beta_2)=sum_(n=-M)^N c_n (E,beta_1) beta_2^n $]
#eq[$ P(beta_2)=beta_2^M f=c_N product_(j=1)^(M+N)(beta_2-beta_(2,j)) $]
#two([
  #label[Scan one angle]
  Fix $E$ and $mu_1$, then let $beta_1=e^(mu_1+i theta_1)$ traverse a circle.
  #v(12pt)
  Solve $P=0$ to obtain the transverse modes at each angle.
],[
  #label[Retain the root geometry]
  #eq[$ ell_j (theta_1)=ln abs(beta_(2,j)(theta_1)) $]
  Moduli encode localization; arguments retain the phases needed to reconstruct the GBZ.
])
#v(1fr)
#source[Generic finite nonzero roots and c_N ≠ 0 shown here; degree loss and multiple roots require separate treatment.]
#speaker-note[1:15. This factorization is elementary but it is the central representation. M is the transverse Laurent pole order; K=M+N counts roots with multiplicity. For matrix Hamiltonians use the determinant polynomial, not the hopping range of a single matrix entry.]

== Amoeba bridge I: integrate one phase analytically

#eq[$ R_E (mu_1,mu_2)=frac(1,(2 pi)^2) integral_0^(2 pi) integral_0^(2 pi) ln abs(f) dif theta_2 dif theta_1 $]
#definition-block([Jensen's identity for each root],[
  $ frac(1,2 pi) integral_0^(2 pi) ln abs(e^(mu_2+i theta_2)-beta_(2,j)) dif theta_2 = max(mu_2,ell_j) $
])
#eq[$ R_E = frac(1,2 pi) integral_0^(2 pi)
  [ln abs(c_N)-M mu_2+sum_j max(mu_2,ell_j (theta_1))] dif theta_1 $]
#v(10pt)
#text(fill:ink)[The two-phase logarithmic integral becomes a one-angle expression in polynomial roots.]
#v(1fr)
#source[Root-factorization/Jensen reduction of the Ronkin definition; all coefficients and roots depend on E and β₁.]
#speaker-note[1:45. All phase integrals in the first line run from zero to two pi. Insert the factorization from the previous slide. The pole gives minus M mu2; each factor gives a maximum by Jensen. Keep the leading-coefficient term: it is essential for the full Ronkin function and its mu1 derivative. The reduction is exact under the displayed generic assumptions, with singular cases treated by limits.]

== Amoeba bridge II: stationarity becomes root counting

#eq[$ a_2=partial_(mu_2) R_E = frac(1,2 pi) integral_0^(2 pi)
  [\# {j:ell_j (theta_1)<mu_2}-M] dif theta_1 $]
#two([
  #label[Only crossings change the count]
  Solve $ell_j (theta_1)=mu_2$. Between crossing angles, the number of enclosed roots stays constant.
  #v(12pt)
  #eq[$ a_2=sum_a (n_a-M) frac(Delta theta_(1,a),2 pi) $]
],[
  #label[Select the Ronkin minimum]
  Search for $a_2=0$, then impose $a_1=0$ through the outer radius search.
  #v(12pt)
  Exchanging axes gives the corresponding root-count identity for $a_1$.
])
#v(1fr)
#source[Derivative formula applies at differentiable points; flat tracks and singular minima need one-sided or subgradient conditions.]
#speaker-note[1:30. Differentiate the maximum away from crossings: each enclosed root contributes one, and the Laurent pole contributes minus M. The angular widths supply exact weights once events are resolved. Exchanging axes explains a1 mathematically; this does not assert that the implementation launches an independent axis-swapped solve. Persistent level equality needs one-sided counts.]

== SGBZ bridge I: transverse OBC becomes a root tie

#eq[$ abs(beta_(2,(1))) <= dots <= abs(beta_(2,(M+N))) $]
#definition-block([The transverse 1D GBZ condition],[
  $ abs(beta_(2,(M)))=abs(beta_(2,(M+1))) $
])
#two([
  #label[Physical balance]
  The two middle modes have equal exponential weight, allowing the transverse OBC bulk-band condition.
],[
  #label[Numerical condition]
  Find zeros of the middle-root gap:
  #eq[$ g(theta_1)=ell_((M+1))(theta_1)-ell_((M))(theta_1) $]
])
#v(1fr)
#source[C. Wang et al., arXiv:2506.22743v2, Sec. II; parenthesized indices denote modulus rank, with multiplicity.]
#speaker-note[1:30. Recall the 1D condition, then apply it at fixed longitudinal beta1. The gap is nonnegative, so a simple sign-change test of the sorted gap is not sufficient. Individual tracked-root differences locate crossings. This transverse condition supplies candidates; the longitudinal strip selection is still required.]

== SGBZ bridge II: replace the strip by loop winding

#eq[$ mu_(2,"mid")(theta_1)=frac(ell_((M))(theta_1)+ell_((M+1))(theta_1),2) $]
#note[Between middle roots, this reference encloses exactly M roots and cancels the transverse Laurent pole.]
#eq[$ w(theta_2)=frac(1,2 pi) integral_0^(2 pi)
 Im[frac(1,f) frac(dif f,dif theta_1)] dif theta_1 $]
#eq[$ W(E,mu_1)=frac(1,2 pi) integral_0^(2 pi) w(theta_2) dif theta_2
 =sum_a w_a frac(Delta theta_(2,a),2 pi) $]
#v(8pt)
#text(fill:ink)[Middle-pair events partition the phase circle into intervals of constant loop winding.]
#v(1fr)
#source[C. Wang et al., arXiv:2506.22743v2, Supplement S2; β₂ = exp[μ₂,mid(θ₁) + iθ₂] along each loop.]
#speaker-note[1:45. The midpoint is a convenient representative inside the middle gap; it is not an additional physical law. Use the equivalent strip-winding formulation in Supplement S2. The derivative is total along the moving reference loop. Each nonsingular loop has integer winding, but the phase average W need not be integer. A reference winding plus event changes determines the interval values.]

== SGBZ selection becomes a scalar radius problem

#definition-block([Longitudinal selection],[
  $ W(E,mu_1^star)=0 $
])
#two([
  #label[Find the localization radius]
  Vary $mu_1$ and evaluate the root-based winding. A bracketed scalar search locates its zero.
  #v(12pt)
  Persistent middle-root ties require one-sided winding limits.
],[
  #label[Recover the Bloch waves]
  At the selected radius, retain the middle-pair roots and their phases.
  #v(12pt)
  A zero-winding plateau still requires a spectral-membership check.
])
#v(1fr)
#source[Project design: winding selection, followed by plateau classification and reconstruction.]
#speaker-note[1:00. We have replaced the large-strip construction by polynomial roots plus a one-dimensional radius search. A vanishing scalar is not automatically a physical spectral point: zero-free plateaus must be distinguished. In a continuum case use the sign of one-sided limits instead of evaluating a singular loop.]

= Numerical design

== Shared roots, different physical constraints

#image("Figures/solver-flow.svg", width:100%)
#v(18pt)
#two([
  #label[Amoeba: constant level]
  Intersections with a trial $mu_2$ determine Ronkin-gradient constraints.
],[
  #label[SGBZ: moving middle gap]
  Middle-root ties and loop winding determine the strip constraint.
])
#v(1fr)
#source[One representation serves both constructions without identifying their physical selection rules.]
#speaker-note[1:00. This architecture follows from the mathematics just derived. The reusable object is a family of transverse roots. The two formulations ask different questions of that family. Avoid discussing class names, interpolation containers, or backend choices.]

== Why follow roots rather than sample a dense grid?

#two([
  #label[Preserve connected modes]
  Root continuation retains branch identity and phase as the longitudinal angle changes.
  #v(14pt)
  Root identity and instantaneous modulus rank are different notions.
],[
  #label[Resolve changes where they occur]
  Adapt the sampling near crossings and degeneracies. Refine the events that change counts or winding.
  #v(14pt)
  Between events, use the constant topological information.
])
#v(18pt)
#definition-block([Design advantage],[Reduce singular quadrature and broad parameter scans to polynomial solves, event locations, and scalar constraints.])
#v(1fr)
#source[This is a structural numerical advantage; quantitative speed and accuracy gains require measured convergence tests.]
#speaker-note[1:15. Root continuation is a means of preserving the geometry, not the main scientific object. Adaptive sampling should be justified by resolving event structure. No timing benchmark is claimed. Finite refinement budgets can still miss close events, so residual and convergence checks remain necessary.]

== Isolated modes, degeneracies, and continuous families

#image("Figures/root-events.svg", width:100%)
#v(12pt)
#note[Equal decay lengths need not mean coincident complex Bloch factors. Preserve the phases of all participating modes.]
#v(1fr)
#source[Conceptual schematic. No computed root trajectory is shown.]
#speaker-note[1:00. Distinct phases can share a modulus. A multiple root is a stronger event and invalidates the regular implicit root derivative. Persistent equality produces line families at fixed E. This is why a numerical design restricted to isolated points is insufficient.]

== Separate the three questions in a GBZ calculation

#table(columns:(1fr,1.75fr),stroke:none,inset:(x:12pt,y:15pt),
  [*Constraint*],[Does the winding or minimum condition hold?],
  [*Spectral membership*],[Does the solution support zeros, or lie in a zero-free plateau?],
  [*Bloch-wave geometry*],[Which phases and connected point or line families survive?],
)
#v(18pt)
#text(fill:ink)[Keep the full complex Bloch factors: the same decay rate can support different spatial oscillations.]
#v(1fr)
#source[Fixed-energy subsets become a sampled GBZ only after an energy sweep.]
#speaker-note[1:00. Tie this slide back to the central hole in the introduction. Solving a scalar equation, deciding spectral membership, and reconstructing a mode family are separate operations. A fixed-energy result is a slice, not an entire GBZ surface.]

= Physical validation and outlook

== An analytic anchor: separable Hatano–Nelson

#eq[$ h=sum_(j=1)^2 J_j (e^(gamma_j) beta_j^(-1)+e^(-gamma_j) beta_j), quad J_j>0 $]
#two([
  #label[Axis-aligned open boundaries]
  #eq[$ beta_j=e^(gamma_j+i k_j) $]
  #eq[$ E=2J_1 cos k_1+2J_2 cos k_2 $]
  A diagonal similarity transform removes the hopping asymmetry.
],[
  #label[Physical checks]
  Recover $mu_j=gamma_j$ and the real spectral interval $[-2(J_1+J_2),2(J_1+J_2)]$.
  #v(12pt)
  Check the reconstructed mode families and their residuals under tighter tolerances.
])
#v(1fr)
#source[Analytic benchmark, real J and γ in the stated basis. No numerical data are plotted.]
#speaker-note[1:15. Substitute beta_j=e^gamma_j z_j. This fixes the sign convention of the skin exponent and gives an exact spectral check. A separable example is an anchor; it does not test generic geometry dependence by itself.]

== What numerical evidence should test the bridge?

#set text(size:20pt)
#table(columns:(1.1fr,1.5fr,1.45fr),stroke:none,inset:(x:10pt,y:13pt),
  table.header([*Physical question*],[*Required data*],[*Comparison*]),
  [Localization],[Model and complex Bloch factors],[Exact benchmark; residuals],
  [Geometry dependence],[Same energy grid, chosen strip bases],[SGBZ versus amoeba support],
  [Mode families],[Roots at fixed E and radius],[Events and connected branches],
  [Numerical reliability],[Repeated tolerance/grid runs],[Stable support; explicit failures],
)
#v(12pt)
#note[Select the Hamiltonian, parameters, lattice basis, diagnostic energy/radius, and energy grid before producing these Figures.]
#v(1fr)
#source[DATA_REQUIREMENTS.md; collect_slide_data.py and plot_slide_data.py produce real-data SVG/PDF Figures.]
#speaker-note[1:00. These are proposed experiments, not results. Use the same model and grid for comparisons. Finite OBC comparisons also require a stated geometry and size sequence. The plotting scripts already exist, and no surrogate or invented dataset is used in this deck.]

== From roots to non-Bloch band geometry

#two([
  #label[Immediate output]
  Spectral support and the complex momenta of the selected mode families.
  #v(16pt)
  Compare decay lengths, phases, and their dependence on strip direction.
],[
  #label[Next physical application]
  Reconstruct connected bands across energies, then study non-Bloch band geometry.
  #v(16pt)
  Eigenvector observables additionally require the Hamiltonian and left/right eigenvectors.
])
#v(1fr)
#source[Band reconstruction remains an experimental extension; the characteristic polynomial alone does not determine eigenvectors.]
#speaker-note[0:45. Motivate the project beyond plotting a spectrum. Avoid presenting a Chern number as a current result. Band separation, completeness and orientation need checks before a topological calculation.]

== The central message

#v(12pt)
#text(size:31pt)[Physical GBZ prescriptions become numerically accessible through polynomial roots.]
#v(24pt)
#two([
  #label[Amoeba]
  Ronkin integral → root moduli → root counts → minimum condition.
],[
  #label[SGBZ]
  Strip boundary problem → middle-root balance → loop winding → radius selection.
])
#v(18pt)
#text(fill:ink)[The project shares the root representation and preserves the distinction between the two physical formulations.]
#speaker-note[0:45. End on the bridge, not implementation details. Invite discussion of which geometry-dependent model most strongly tests the two selections.]

== References

#set text(size:20pt)
H.-Y. Wang, F. Song, and Z. Wang, *Amoeba Formulation of Non-Bloch Band Theory in Arbitrary Dimensions*, Phys. Rev. X *14*, 021011 (2024).\
#link("https://doi.org/10.1103/PhysRevX.14.021011")
#v(24pt)
C. Wang et al., *Universal theory for geometry-dependent non-Hermitian bands*, arXiv:2506.22743v2 (2025).\
#link("https://arxiv.org/abs/2506.22743v2")
#v(24pt)
#note[Derivation assumptions, source mapping, and numerical data requirements accompany the editable deck.]
#speaker-note[References for discussion. The slide equations use the project's beta=exp(mu+i theta) convention. See sources.md for theory versus numerical-design attribution.]

// Optional computed Figures, only after real data have been supplied.
#if sys.inputs.at("with-data", default: "false") == "true" [
= Computed data
== Computed root tracks
#align(center)[#image("Figures/root-tracks.svg", height: 365pt)]
#source[Model and calculation provenance: Figures/plot-provenance.json]
== Computed subsets on the angle torus
#align(center)[#image("Figures/subsets.svg", height: 365pt)]
#source[Model and calculation provenance: Figures/plot-provenance.json]
== Computed energy sweep
#align(center)[#image("Figures/spectrum.svg", height: 365pt)]
#source[Model and calculation provenance: Figures/plot-provenance.json]
]

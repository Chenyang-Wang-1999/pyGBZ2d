#import "style.typ": *  // modified from @preview/zhaji:0.1.0

#let nt = note(
  title: "Chern number calculation: QWZ benchmark",
  author: "Chenyang Wang",
  font-head: "Arial",
	font-text: "Times New Roman",
	first-line-indent: 0em,
	lang: "en"
)

#show: nt.make
#show raw.where(block: true): it => block(it, fill: rgb("edf6fd"), inset: 5pt, breakable: true)
#show raw: set text(size: 8.5pt)
#show raw.where(block: false): it => highlight(it, fill: rgb("e0e0e0"), top-edge: 1.2em, bottom-edge: -.5em, extent: 0.2em)

#let mathbf(x) = $bold(upright(#x))$
#let rme = $upright(e)$
#let rmi = $upright(i)$

= Chern number calculation

#let results = json("Figures/Chern-number-calculation/benchmark-summary.json")
#let fmt(x, digits: 6) = str(calc.round(x, digits: digits))
#let fig(name) = image("Figures/Chern-number-calculation/" + name + ".svg", width: 100%)
#let chern-rows = results.at("chern", default: ())

This benchmark applies the `pyGBZ2d` solver and experimental mesh/Chern pipeline to the standard non-Hermitian Qi-Wu-Zhang (QWZ) model @yao2018nonhermitianchern @wang2024amoeba. Right-state integration and the corrected biorthogonal integration both reproduce the literature's lower-band Chern numbers: *one for the nontrivial case and zero for the trivial case*. The biorthogonal correction removes a second-order link error on the same saved meshes. Its derivation and numerical checks are given below; a few mesh edges still exceed the refinement target.

The #link("Chern-number-calculation.py")[companion script] runs a $20 times 20$ coarse energy sweep, a $201 times 201$ fine sweep, adaptive mesh refinement, and Chern integration. Data default to `application/data/QWZ-benchmark`; sweeps and refinement solves use 20 worker processes.

== Model and literature targets

In units $t=v=1$, the Bloch Hamiltonian is
$ h(k_x,k_y)=(sin k_x+rmi gamma)sigma_x+(sin k_y+rmi gamma)sigma_y
 +(m-cos k_x-cos k_y)sigma_z. $
For complex Bloch factors, use
$ sin k_j=(beta_j-beta_j^(-1))/(2rmi), quad
cos k_j=(beta_j+beta_j^(-1))/2, quad beta_j=exp(mu_j+rmi theta_j). $
The logarithmic radii $mu_j$ are solved numerically. Taking only `angle(beta)` would discard the skin decay and evaluate a different Hamiltonian.

#table(columns: (1fr, auto, auto, 1fr), inset: 6pt, stroke: .4pt + luma(75%),
  table.header([Case], [$m$], [$gamma$], [Literature lower-band $C$]),
  [Nontrivial], [1.2], [0.4], [1],
  [Trivial], [2.6], [0.4], [0],
)

The parameters and targets follow Fig. 9 of Ref. @wang2024amoeba. The lower band has $"Re" E<0$; sign conventions are reconciled below. The reported phase boundary $m_c=2sqrt(1+gamma^2)$ supplies context but is not used to assign computed Chern numbers.

=== Characteristic polynomial

BerryPy builds the two-orbital square-lattice model and generates its Laurent polynomial. For these parameters it also has a compact independent expression:
$ f(E,beta_x,beta_y)= & E^2-(m^2+2-2gamma^2) \
 & +(m-gamma)(beta_x+beta_y)+(m+gamma)(beta_x^(-1)+beta_y^(-1)) \
 & -1/2(beta_x+beta_x^(-1))(beta_y+beta_y^(-1)). $
The script checks this ten-term polynomial against $det(E I-h)$ at eight deterministic complex arguments per case, then calls `pygbz2d.amoeba.collect_GBZ_subsets` with default numerical settings and the NumPy backend.

#pagebreak()
== Energy-window discovery and the fine sweep

The coarse real-energy range is $[-m-2.2,m+2.2]$. The 20-point imaginary grid includes the real axis and encloses $[-sqrt(2)gamma,sqrt(2)gamma]$. Including zero matters because a thin spectrum may be missed by an even grid symmetric about, but not containing, zero.

The coarse calculation found 36 spectral energies for $m=1.2$ and 16 for $m=2.6$, with no failed nodes. These observations locate the main spectrum but do not certify its boundary. In fact, an initial narrow fine window still contained spectral points on its outer row: the coarse grid missed a thin complex-energy lobe. That partial calculation was retained separately and excluded from the benchmark mesh.

The corrected fine window retains two coarse spacings of real-energy margin and a conservative imaginary enclosure. Since all intercell hoppings are Hermitian and $(H-H^dagger)/(2rmi)=gamma(sigma_x+sigma_y)$ on every site, every finite OBC eigenvalue obeys $abs("Im" E)<=sqrt(2)gamma$. The scan extends one coarse imaginary spacing beyond this bound. The final grid boundary is checked again before triangulation.

#figure(fig("sweeps"), caption: [
  Coarse and fine energy scans. Colored points are successful GBZ results with nonempty subsets; gray points are successfully classified exterior points. Failed solves, if any, are marked separately. Completed fine panels zoom in vertically on the spectrum; the table gives the full computed windows. Coarse data are used only for window discovery; the GBZ meshes use the fine data.
])

#table(columns: (auto, 1fr, 1fr, auto), inset: 5pt, stroke: .4pt + luma(75%),
  table.header([$m$], [Fine Re interval], [Fine Im interval], [Completed]),
  ..("nontrivial", "trivial").map(key => {
    let r=results.at(key).at("fine", default: none)
    if r == none { ([#results.at(key).mass], [Pending], [Pending], [0 / 40401]) }
    else { ([#results.at(key).mass], [#fmt(r.real_range.at(0), digits: 3) to #fmt(r.real_range.at(1), digits: 3)],
      [#fmt(r.imag_range.at(0), digits: 3) to #fmt(r.imag_range.at(1), digits: 3)], [#r.done / #r.total]) }
  }).flatten(),
)

#pagebreak()
== From GBZ subsets to an integration surface

The solver returns `GBZResult` objects. Every `PointSubset` contributes one sample and every stored row of a `LineSubset` contributes one sample. The application uses `flatten_results` to collect $(E,mu_x,theta_x,mu_y,theta_y)$ and retains the energy-slice provenance. A failed solve is never relabeled as a spectral exterior point.

For this two-band benchmark, the real line gap permits a direct split by the sign of $"Re" E$. This avoids selecting a radius-graph clustering threshold merely to identify the two already separated bands. Within each band, all samples are passed to the library's `build_band_mesh` routine. It merges coincident phase vertices and constructs a periodic Delaunay triangulation on the phase torus.

The next stage calls `pygbz2d.experimental.refine_mesh`. It builds an energy interpolant to predict energies at long-edge midpoints, calls the actual GBZ solver at those energies, and accepts matched computed GBZ samples. Index-transition fallback probes are handled by the existing refinement implementation. Predicted values alone are not inserted as solved eigenstates. The requested maximum phase-space edge length is $0.2$ radians, with at most ten outer refinement rounds; each batch uses 20 workers.

For this trace-zero two-band model, $f(E,beta)=f(-E,beta)$, so the Ronkin minimizers and GBZ coincide for partner energies. The lower-band torus is refined once per parameter set and reused for the upper band with $E mapsto -E$. Both bands' eigenvectors and Chern integrals are evaluated independently. This reduces duplicate geometric work; all nodes of the original energy sweeps are still solved directly.

The report compares the initial fine-sweep mesh with the refined mesh. It records vertex and triangle counts, maximum and 99th-percentile edge lengths, connected components, Euler characteristic, and the incidence of every edge. For a single closed torus one expects one component, $chi=V-E+F=0$, and exactly two incident faces per edge. These combinatorial tests must be considered together with sampling density and eigenvector residuals: triangulating a sparse cloud can formally close a hole without having resolved the underlying GBZ there.

=== Chern-number convention

At every mesh vertex the routine computes null vectors of $E I-h(beta_x,beta_y)$. Right vectors are normalized to unit Euclidean norm; in biorthogonal mode left vectors additionally obey $⟨u_L|u_R⟩=1$. Right-state (RR) flux remains $-arg(⟨R_i|R_j⟩⟨R_j|R_k⟩⟨R_k|R_i⟩)$. Biorthogonal (LR) flux now uses the reciprocal transports derived on the following pages. Triangles are oriented positively in $(theta_x,theta_y)$.

The source's raw-sign result is $C_"raw"=sum "Re" Phi/(2pi)$. The literature convention used here is
$ C_"lit"=-C_"raw", $
as fixed by the known Hermitian QWZ limit and the sign of the Berry connection in Ref. @wang2024amoeba. Both values are saved; this conversion applies uniformly to all bands and modes. Here “raw” labels the sign convention, not the obsolete raw-overlap LR rule. Eigenvector residuals, nullities, triangle fluxes, and link validity are checked separately.

#pagebreak()
== Biorthogonal curvature and the local link error

Consider one smooth, nondegenerate band over real coordinates $q=(q_1,q_2)$ on the integration surface, with $⟨L|R⟩=1$. On the GBZ these derivatives include the variation of $mu(theta)$. Define
$ a_alpha=⟨L|partial_alpha R⟩, quad A_alpha=rmi a_alpha, quad P=|R⟩⟨L|, $
$ B_(1 2)=partial_1 A_2-partial_2 A_1
 =rmi(⟨partial_1 L|partial_2 R⟩-⟨partial_2 L|partial_1 R⟩). $
The antisymmetrization exchanges the coordinate indices; both terms retain a left bra and a right ket. Shen, Zhen, and Fu establish equality of the integrated RR, LL, LR, and RL Chern numbers @shen2018nonhermitiantopology. Their local curvatures can differ. The normalization conditions and these four definitions are also reviewed by Ashida, Gong, and Ueda @ashida2020nonhermitian.

For a sufficiently small oriented loop,
$ W=exp(integral.cont a)=exp(-rmi integral B dif q_1 dif q_2), quad
Phi=rmi "Log" W=-arg W+rmi log abs(W). $
Thus the Wilson-loop description remains valid for biorthogonal states, including a generally complex local curvature. What needs correction is the finite-edge approximation to parallel transport.

=== Second-order expansion of the overlap

Let $U_(i j)=⟨L(q_i)|R(q_j)⟩$ and define
$ Q_(alpha beta)=⟨partial_alpha L|(1-P)|partial_beta R⟩, quad
g_(alpha beta)=(Q_(alpha beta)+Q_(beta alpha))/2. $
For a straight edge of displacement $delta$, Taylor expansion gives, with summed coordinate indices,
$ "Log" U(q,q+delta)
 =a_alpha delta_alpha+1/2(partial_alpha a_beta-Q_(alpha beta))delta_alpha delta_beta+O(abs(delta)^3) $
$ "Log" U(q,q+delta)
 =integral_q^(q+delta) a-1/2 g_(alpha beta)delta_alpha delta_beta+O(abs(delta)^3). $
This follows from $partial_alpha a_beta=⟨partial_alpha L|partial_beta R⟩+⟨L|partial_alpha partial_beta R⟩$ and $⟨partial_alpha L|R⟩=-a_alpha$.

Around a triangle of diameter $h$ and area of order $h^2$, the old phase therefore contains
$ -arg product_(partial triangle) U
 ="Re" integral_triangle B+1/2 "Im" sum_(partial triangle)
 g_(alpha beta)delta_alpha delta_beta+O(h^3). $
The symmetric tensor $g$ is real in the Hermitian case, but generally complex for LR states. Its imaginary part contaminates the local phase at the same order as the physical flux. Dividing by area does not remove this error, and merely replacing each $U$ by $U/abs(U)$ leaves the erroneous phase intact.

Raw overlaps can converge along a *fixed contour* whose number of subdivisions increases. Shrinking a triangle while retaining only three links is a different limit. The latter needs an edge approximation accurate beyond second order.

#pagebreak()
== Reciprocal biorthogonal links

For each undirected edge, evaluate both LR overlaps and set
$ s_(i j)=sqrt(U_(i j)U_(j i)), quad
T_(i j)=U_(i j)/s_(i j), quad T_(j i)=T_(i j)^(-1). $
Both directions use the *same square root*, chosen continuously from $+1$ for short edges. Since
$ "Log"(U_(i j)U_(j i))=-g_(alpha beta)delta_alpha delta_beta+O(h^3), $
the symmetric second-order term cancels, giving
$ "Log" T_(i j)=integral_i^j a+O(h^3), quad
Phi_(i j k)=rmi "Log"(T_(i j)T_(j k)T_(k i)). $
For shape-regular triangles and bounded derivatives this gives an $O(h^3)$ flux error, or $O(h)$ curvature error after division by area. In the Hermitian limit, $U_(j i)=U_(i j)^*$ and $T_(i j)=U_(i j)/abs(U_(i j))$.

Under $|R_i⟩ mapsto c_i|R_i⟩$ and $⟨L_i| mapsto c_i^(-1)⟨L_i|$, for arbitrary nonzero complex $c_i$, the product $U_(i j)U_(j i)$ is invariant and $T_(i j) mapsto (c_j/c_i)T_(i j)$. Loop products are invariant. Independently taking the principal square root of $U_(i j)/U_(j i)$ would introduce gauge-dependent signs and is not the implemented prescription.

Fukui, Hatsugai, and Suzuki use normalized Hermitian links and explicit *inverse links for backward traversal* in their lattice Chern construction @fukui2005chern. Their work establishes gauge invariance, lattice integer quantization, and the role of admissibility. The LR square-root normalization above is derived here; it is not attributed as an equation from that paper.

=== Implementation and validity checks

`experimental.chern` evaluates each canonical undirected edge once and sums its signed phase and log magnitude around each triangle. This enforces inverse backward transport and avoids multiplying potentially large or small loop factors. The saved `flux` remains real; `flux_complex` retains both parts, and `link_rule` identifies `biorthogonal-reciprocal-v1`.

A zero or nonfinite left/right self-overlap now raises an error rather than skipping normalization. Zero/nonfinite edge products and products on the negative-real square-root branch cut are rejected. The principal root is used away from the cut, as the continuation of the root near one for sufficiently short edges. These checks do not certify band smoothness on a coarse edge.

On a closed consistently oriented mesh, reciprocal edges cancel pairwise: the total imaginary flux vanishes, and the summed principal real flux is an integer multiple of $2pi$ up to roundoff. Recovering the continuum integer also requires resolved geometry and triangle phases. An integer alone is not a convergence certificate.

=== Analytic regression examples

For $|R⟩=(1,x)^T$ and $⟨L|=(1-s x,s)$ with $s=1+rmi$, the smooth Hamiltonian $H=2P-I$ has eigenvalues $plus.minus 1$ and zero curvature. On vertices $(0,0),(h,0),(0,h)$ the old overlap loop is $1-2rmi h^2$, giving an apparent curvature tending to four. The corrected loop gives zero.

Replacing $s$ by $(1+rmi)y$ gives $a=(1+rmi)y dif x$ and $B=1-rmi$. The corrected flux divided by triangle area converges to this full complex value. Tests also cover irregular QWZ tori in both phases, complex gauge transformations, vertex relabeling, orientation reversal, the Hermitian limit, invalid links, and saved complex flux. The full derivation and API notes are in `doc/experimental.md`.

#pagebreak()
== Numerical results and checks

Both fine sweeps finished with zero failed nodes and zero spectral hits on the outer grid. They contain 4,302 and 1,416 spectral energies, respectively, producing 25,836 and 8,440 GBZ samples. Maximum normalized polynomial residuals are $6.23 times 10^(-14)$ and $6.57 times 10^(-15)$.

=== Independent BZ and Hermitian checks

Before the GBZ calculation, the same mesh/Chern routines were evaluated on regular $21 times 21$, $41 times 41$, and $81 times 81$ phase grids, both at $gamma=0$ and at $gamma=0.4$. For $m=1.2$ the lower-band result is $C_"lit"=1$; for $m=2.6$ it is zero. Both right-state and biorthogonal modes agree to floating-point precision, with no near-zero links. These are checks of modeling, integration, and sign conventions; they do not replace the GBZ calculation.

=== Chern numbers from the computed GBZ

Completed adaptive-refinement outputs: #results.at("refined_band_meshes", default: 0) / 4 band meshes. This count records returned meshes; the edge-length and topology checks determine whether they meet the requested quality.

#if chern-rows.len() == 0 [
  The fine energy sweeps are in progress. No GBZ Chern number is reported yet.
] else {
  table(columns: (auto, auto, 1fr, auto, auto), inset: 5pt, stroke: .4pt + luma(75%),
    table.header([$m$], [Band], [Mesh / mode], [$C_"lit"$], [Max. edge]),
    ..chern-rows.map(r => ([#r.mass], [#r.band],
      [#if r.level == "before refinement" [Base] else [Refined] / #if r.mode == "right" [RR] else [LR]],
      [#fmt(r.literature_sign_chern, digits: 8)], [#fmt(r.max_edge, digits: 3)])).flatten())
}

RR denotes right-state mode and LR denotes biorthogonal mode. Values are displayed to eight decimal places; `benchmark-summary.json` retains unrounded results, raw-sign values, and full numerical diagnostics. The literature target for the lower band is $1$ at $m=1.2$ and $0$ at $m=2.6$; a benchmark conclusion requires both agreement with that target and satisfactory mesh/eigenvector checks.

#if chern-rows.any(r => r.level == "after refinement") {
  pagebreak()
  heading(level: 2)[GBZ meshes and Berry flux]
  figure(fig("gbz-and-flux"), caption: [
    Lower-band meshes after adaptive refinement, colored by $mu_x$, and their right-state triangle Berry fluxes in the literature sign convention. Red edges exceed the 0.2-radian target. Seam-crossing triangles are hidden only in this flat visualization; they remain present in the Chern integral. The lower panels show flux per triangle, rather than curvature per unit area.
  ])
  table(columns: (auto, auto, 1fr, 1fr, auto, 1fr), inset: 6pt, stroke: .4pt + luma(75%),
    table.header([$m$], [Band], [Vertices], [Triangles], [$chi$], [Max. edge]),
    ..chern-rows.filter(r => r.level == "after refinement" and r.mode == "right").map(r => (
      [#r.mass], [#r.band], [#r.n_vertices], [#r.n_triangles], [#r.topology.chi], [#fmt(r.max_edge, digits: 6)]
    )).flatten())
  [
    Partner bands share the refined geometric surface. The full JSON records edge incidence, unused vertices, connected components, and total phase area, as well as the eigenvector diagnostics. A closed torus and an integer right-state result do not, by themselves, certify the requested maximum edge length.
  ]
}

#pagebreak()
== Mesh limits and the resolved integration defect

=== Remaining mesh-quality limitations

For $m=1.2$, ten refinement rounds increase the lower-band mesh from 12,918 to 15,274 vertices. Eight unique edges near the thin complex-energy tips still exceed 0.2 radians; the longest is 0.23192382 radians. The global edge target has not been reached.

The $m=2.6$ mesh grows from 4,220 to 5,841 vertices but reaches the ten-round cap with 23 long edges and maximum edge 0.25348215 radians. Both meshes remain closed, consistently oriented tori of area $4pi^2$, with no unused vertices or nonmanifold edges and one-dimensional eigenvector null spaces throughout.

The final rounds increasingly merge nearby solved points into existing vertices, leaving the longest edge unchanged; target-match distances reach about 0.188 radians. These refinement limits remain separate from the link correction. An integer Chern number does not certify geometric convergence.

=== Before and after the LR correction

Before correction, the nontrivial refined mesh returned $C_"lit"=1$ in RR mode but approximately $1.00589930$ with the raw LR overlap triangle rule. The saved historical values show its orientation defect:

#table(columns: (auto, 1fr, 1fr, 1fr), inset: 6pt, stroke: .4pt + luma(75%),
  table.header([$m$], [Historical mode], [Positive orientation], [Negative orientation]),
  [1.2], [Right-state], [$-1.00000000$], [$+1.00000000$],
  [1.2], [Biorthogonal], [$-1.00589930$], [$+0.99410070$],
  [2.6], [Right-state], [$approx 0$], [$approx 0$],
  [2.6], [Biorthogonal], [$+1.87413 times 10^(-5)$], [$+1.87413 times 10^(-5)$],
)

For the raw LR overlaps, the product of all triangle loops retains factors $U_(i j)U_(j i)$ from internal edges. Their measured contribution in the source sign convention is
$ B_"edge"=-1/(2pi) sum_({i,j}) arg(U_(i j)U_(j i))=-0.005899300409245. $
This accounts for the old fractional offset and the even part under orientation reversal. For $m=2.6$ the contribution is $+1.87413 times 10^(-5)$. Nearly symmetric base meshes cancel it almost exactly, explaining why regular-grid tests alone did not expose the defect. Eigenvector residuals remain below $2.3 times 10^(-13)$.

Reintegration with reciprocal links uses the same eight input mesh files, verified by SHA-256. The updated main results table uses this corrected rule. The lower-band orientation checks now give:

#table(columns: (auto, 1fr, 1fr, 1fr), inset: 5pt, stroke: .4pt + luma(75%),
  table.header([$m$], [Corrected mode], [Positive orientation], [Negative orientation]),
  ..("nontrivial", "trivial").map(key => {
    let checks=results.diagnostics.at(key).orientation_checks.filter(r => r.mode == "biorthogonal")
    ([#results.at(key).mass], [LR reciprocal], [#fmt(checks.at(0).raw_chern, digits: 10)],
     [#fmt(checks.at(1).raw_chern, digits: 10)])
  }).flatten(),
)

The maximum per-triangle complex gauge error is $5.84 times 10^(-15)$ over both refined lower-band meshes. Reversing triangle order gives exactly opposite complex fluxes in the local check; independent orientation runs agree up to roundoff. The minimum real parts of the invariant edge products are 0.98694 and 0.97094, safely away from the square-root cut. The imaginary parts of the integrated lower-band Chern numbers are below $3 times 10^(-18)$. `diagnose` reproduces these checks without new GBZ solves. The earlier outputs remain in `biorthogonal-raw-overlap-baseline`.

#pagebreak()
== Reproduction and data files

Use an environment with NumPy, SciPy, Matplotlib, BerryPy, and SymPy, and install this checkout with `python -m pip install -e .`. Run from the repository root:

```sh
python application/Chern-number-calculation.py coarse --workers 20
python application/Chern-number-calculation.py fine --workers 20
python application/Chern-number-calculation.py mesh --workers 20
python application/Chern-number-calculation.py chern
python application/Chern-number-calculation.py bz-check
python application/Chern-number-calculation.py diagnose
python application/Chern-number-calculation.py report
typst compile application/Chern-number-calculation.typ
```

`all` runs coarse, fine, and mesh stages. Both cases are selected by default; `--cases` selects one. `--data-dir` defaults to `application/data/QWZ-benchmark`. Compatible sweep checkpoints resume; incompatible grids or model parameters are rejected. Every energy-grid node is solved directly.

The optional `chern` stage reintegrates saved meshes without rewriting them or rerunning scans/refinement. It was used for this correction. `mesh` already integrates its outputs. Use a fresh directory for changed solver/refinement settings. Partner-band mesh reuse requires this model's even-in-$E$, trace-zero symmetry.

Files distinguish `coarse`, `fine`, `base-mesh`, and refined `mesh` data. Sweeps retain full solver results, polynomial/grid data, and timings. Meshes retain phases, connectivity, energies, and logarithmic radii. JSON retains both Chern conventions, modes, and integration rules. Flux archives contain real/complex fluxes and oriented triangles. The summary includes current source hashes; old outputs and their source version are archived separately.

The partial scan `QWZ-nontrivial-fine-truncated-window.pkl` is retained and excluded from subsequent stages. Sweep errors remain in result objects. Refinement logs report accepted/rejected counts and fallback totals, but do not retain every probe's full result. Solver tolerances are unchanged.

`pygbz2d.amoeba` supplies GBZ results; `pygbz2d.experimental` supplies triangulation, refinement, and integration. The Python example owns file I/O and plotting, uses installed package APIs, and can be copied outside this checkout. No sibling scripts are needed. The dependency audit is in `doc/experimental.md`. Regression tests: `python -m pytest tests/test_chern_biorthogonal.py -q`.

== References
#bibliography("pyGBZ2d.bib", title: none, style: "american-physics-society")


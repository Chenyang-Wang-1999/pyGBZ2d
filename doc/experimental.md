# Experimental algorithms

`pygbz2d.experimental` contains trial-stage algorithms with limited validation.
APIs may change without notice. Numerical robustness and convergence across
all models and parameter regimes have not been established. Passing the
reported benchmarks does not establish general reliability; known limitations
and the checks needed to interpret results are documented below.

These routines are reusable library code, with explicit opt-in imports. They
cover band clustering, periodic band meshes, adaptive refinement and Chern
integration, using NumPy/SciPy and the rest of `pygbz2d`. Package placement
expresses reuse, while the `experimental` namespace records their maturity.

- [Band clustering](#band-clustering)
- [Periodic band meshes](#periodic-band-meshes)
- [Adaptive refinement](#adaptive-refinement)
- [Chern integration](#chern-integration)
- [Parameters and constants](#parameters-and-constants)
- [Application boundary and dependency audit](#application-boundary-and-dependency-audit)
- [Validation and reproduction](#validation-and-reproduction)
- [References](#references)

## Band clustering

`experimental.band_clustering` groups GBZ sweep samples by connected components
of a radius graph. `flatten_results` retains every `PointSubset` and every
stored sample of each `LineSubset`, with energy-slice and subset provenance.
Failed/empty results are skipped, and the `BandPoints.n_dropped` diagnostic
records nonfinite samples. Line samples are not decimated: nearly vertical
arcs in the phase torus need their native sampling to preserve connectivity.

The embedding used by `embed` is

```text
[alpha_E * Re(E)/d_re, alpha_E * Im(E)/d_im,
 w_mu * mu1, w_mu * mu2,
 cos(theta1), sin(theta1), cos(theta2), sin(theta2)]
```

The angular pairs handle the periodic seam. The energy block helps separate
bands whose phase/radius projections overlap, while the radius block retains
separation between states with different decay lengths. By default,
`w_mu = 1/max(ptp(mu1), ptp(mu2), 1)`; omitted energy steps are inferred from
median spacings. Supply the original `d_re` and `d_im` for enriched or
irregular sweeps, because extra probe energies can distort inferred spacings.

```python
from pygbz2d.experimental import cluster_bands, summarize_clusters

cl = cluster_bands(results, d_re=d_re, d_im=d_im)
print(cl.eps, cl.eps_window, cl.n_clusters)
print(summarize_clusters(cl.points, cl.labels))
```

`cluster_bands` accepts exactly one of `results` or a pre-flattened `points`
argument. An explicit `eps` selects the graph radius. Otherwise the algorithm
scans radii and selects the geometric midpoint of the widest eligible
plateau in log-radius. Eligibility requires enough clusters and a largest
component containing a sufficient fraction of the cloud. If no eligible
plateau exists, it uses the first full-merge radius, or the last scanned
radius if the cloud never fully merges. The scan stops at full merging and
has an absolute radius cap to bound pair-graph growth.

Inspect `scan`, `eps_window`, nearest-neighbor percentiles and inter-cluster
`margins` before identifying components with physical bands. A graph component
can split a poorly sampled band or join nearby bands; a stable component
count is evidence of a useful sampling scale, not a proof of band separation.
For a benchmark with an independently known real line gap, the application
may select bands directly by the sign of Re(E), as in the QWZ example below.

| API | Role |
| --- | --- |
| `BandPoints`, `flatten_results` | Samples and source-result provenance |
| `grid_steps`, `embed` | Energy scales and periodic embedding |
| `radius_graph_labels`, `cluster_bands` | Radius-graph connected components |
| `knn_distance_stats`, `eps_stability_scan`, `widest_plateau` | Sampling-scale diagnostics and radius selection |
| `inter_cluster_margins`, `summarize_clusters` | Separation and per-component summaries |
| `BandClustering` | Labels, selected radius, plateau, scan and diagnostics |

## Periodic band meshes

```python
import numpy as np
from pygbz2d.experimental import flatten_results, build_band_mesh

# Select one smooth band before constructing its projected phase surface.
selected = [r for r in results if r.success and r.is_gbz and r.E_ref.real < 0]
points = flatten_results(selected)
verts, triangles, data, n_merged = build_band_mesh(points)
```

`verts` contains `(theta1, theta2)` in radians; `triangles` holds vertex ids.
`data` contains `E`, `mu1`, `mu2` and `slice_idx`, retaining provenance into
`selected`. The optional `mask` selects part of a `BandPoints` cloud. For a
`BandClustering` result, `build_cluster_mesh(cl, cluster_id)` is a thin adapter.

Near-coincident projections are merged, retaining the first sample. The
periodic Delaunay algorithm replicates the phase cloud into neighboring
tiles, uses a translation-consistent tie break, and maps faces back to their
original vertex ids. The original algorithm and numerical tolerances are
retained. It assumes a band can be represented as one sheet over the phase
torus. A closed triangulation of a sparse projection does not prove that the
GBZ is fully resolved or single-valued there.

`mesh_topology(triangles)` reports Euler characteristic, components and
boundary loops, using the vertex-id range `0..max(triangles)`.
`edge_lengths` returns directed face edges with duplicates and their short-arc
lengths; `triangle_areas` uses local periodic unwrapping. The application also
checks unused vertices, edge incidence and total phase area.

## Adaptive refinement

```python
from pygbz2d.experimental import refine_mesh, edge_lengths

indices = np.asarray([selected[s].index for s in data["slice_idx"]], dtype=int)
verts, triangles, data = refine_mesh(
    verts, triangles, data, coeffs, degs, "amoeba", indices,
    edge_thresh=0.2, max_iters=10, n_procs=20, verbose=True,
)
maximum_edge = edge_lengths(verts, triangles)[1].max()
```

`indices` has shape `(n_vertices, 2)` and contains the source results'
`(n_0D, n_1D)` labels. The two-level refinement is unchanged: a periodic
Clough-Tocher interpolant predicts midpoint energies from real samples;
dyadic geometry generates pending positions; a batch of actual GBZ solves
and index-boundary fallback probes supplies accepted eigenstates. Predicted
energies alone never become real mesh data. Returned `data` includes `E`,
`mu1`, `mu2` and, after accepted refinement, `idx0`/`idx1` labels.

The routine returns the current mesh when the edge target is met, no new
solved samples can be accepted, or the iteration cap is reached. Check the
returned edges before claiming convergence. The pre-existing nearest-match
refinement limits observed in the QWZ example are not changed by extraction.
`verbose` enables progress output. The batch worker is defined in the
installed module so Windows spawn workers do not require a demo script.

Omitted `edge_thresh`, `match_tol` and `max_iters` resolve from the module's
live constants. `solve_batch` resolves its bisection budget before creating
workers. See [Parameters and constants](#parameters-and-constants) for defaults.

## Chern integration

```python
from pygbz2d.experimental import integrate_chern

result, flux = integrate_chern(
    verts, triangles, data["E"], data["mu1"], data["mu2"], Hfun,
    mode="biorthogonal", orientation="positive",
)
print(result.chern, result.link_rule)
# Caller-owned persistence, when desired:
np.savez("band-flux.npz", **flux)
```

`Hfun(beta1, beta2)` returns the dense Hamiltonian. `ChernResult` is a scalar
diagnostic dataclass, with nested `VertexDiagnostics`; it contains no file
path. The returned `flux` dictionary contains `flux`, `flux_complex`,
`link_rule`, `oriented_triangles`, `beta1`, and `beta2`.
The LR derivation, reference citations, invalid-link handling and sign
convention follow below.

`pygbz2d.experimental.chern` is an experimental, single-band postprocessor.
It reconstructs `beta_j = exp(mu_j + i theta_j)` from torus-mesh arrays,
extracts eigenvectors of `H(beta1, beta2)` at the stored energy, and integrates
triangle Berry fluxes. The QWZ application and its numerical results are in
[`application/Chern-number-calculation.typ`](../application/Chern-number-calculation.typ).
The numerical API has no file/plotting or playground dependencies; the
application owns file I/O. See [Application boundary and dependency audit](#application-boundary-and-dependency-audit)
for the dependency audit.

### Continuous definition and sign

Assume a smooth, nondegenerate band on an oriented closed two-dimensional
surface, with local frames normalized by $\langle L|R\rangle=1$. Derivatives
below refer to real surface coordinates $q=(q_1,q_2)$; on the GBZ they are
derivatives of the pulled-back states, including the variation of both
$\theta$ and $\mu(\theta)$.

\[
a_\alpha=\langle L|\partial_\alpha R\rangle,\qquad
A_\alpha=i a_\alpha,\qquad P=|R\rangle\langle L|,
\]
\[
B_{12}=\partial_1 A_2-\partial_2 A_1
=i\bigl(\langle\partial_1 L|\partial_2 R\rangle
-\langle\partial_2 L|\partial_1 R\rangle\bigr).
\]

Both terms retain the left bra and right ket; antisymmetrization exchanges
the coordinate indices. For a small oriented loop,

\[
W=\exp\oint a=\exp\left(-i\int B\,dq_1dq_2\right),\qquad
\Phi=i\operatorname{Log} W=-\operatorname{Arg}W+i\log|W|.
\]

The real flux gives `chern = sum(Re Phi)/(2*pi)`. This is the source's raw
sign convention. The QWZ application also records `literature_sign_chern =
-chern`, using the opposite Berry-connection sign in its reference comparison.
Changing the link discretization does not change this sign conversion.

Shen, Zhen, and Fu prove that the four left/right definitions have equal
integrated Chern numbers, despite different local curvatures [1]. The review
by Ashida, Gong, and Ueda discusses these definitions and their normalization
conditions [2]. Equality concerns the continuum invariant on the same smooth
band bundle; it does not validate every finite-mesh overlap formula.

### Why the raw LR triangle fails

Set $U_{ij}=\langle L(q_i)|R(q_j)\rangle$. This is a first-order approximation
to parallel transport, but its second-order term matters for local curvature.
Define

\[
Q_{\alpha\beta}=\langle\partial_\alpha L|(1-P)|\partial_\beta R\rangle,
\qquad g_{\alpha\beta}=(Q_{\alpha\beta}+Q_{\beta\alpha})/2.
\]

From $\partial_\alpha a_\beta=\langle\partial_\alpha L|\partial_\beta R\rangle
+\langle L|\partial_\alpha\partial_\beta R\rangle$ and
$\langle\partial_\alpha L|R\rangle=-a_\alpha$, Taylor expansion on a straight
edge of displacement $\delta$ gives (summed coordinate indices)

\[
\begin{aligned}
\operatorname{Log}U(q,q+\delta)
&=a_\alpha\delta_\alpha+
\tfrac12(\partial_\alpha a_\beta-Q_{\alpha\beta})\delta_\alpha\delta_\beta
+O(|\delta|^3)\\
&=\int_{q}^{q+\delta}a
-\tfrac12g_{\alpha\beta}\delta_\alpha\delta_\beta+O(|\delta|^3).
\end{aligned}
\]

Consequently, on a triangle of diameter $h$ and area of order $h^2$,

\[
-\operatorname{Arg}\prod_{\partial\triangle}U
=\operatorname{Re}\int_\triangle B
+\tfrac12\operatorname{Im}\sum_{\partial\triangle}
g_{\alpha\beta}\delta_\alpha\delta_\beta+O(h^3).
\]

The symmetric tensor $g$ is real for normalized Hermitian states, but it is
generally complex in the LR case. Its phase error is of the same order as
the triangle's flux, so dividing by area does not remove it. Likewise,
$U_{ij}U_{ji}\ne 1$ and its phase need not vanish: shared internal edges
fail to cancel. Normalizing each raw LR overlap by its own absolute value
preserves this erroneous phase.

This does **not** rule out biorthogonal Wilson loops. Raw overlaps converge
along a **fixed contour** when the number of subdivisions increases: the
sum of their second-order errors goes to zero. Shrinking a triangle with
only three edges is a different limit.

### Reciprocal transport used by the implementation

For each undirected edge, evaluate both overlaps and define

\[
s_{ij}=\sqrt{U_{ij}U_{ji}},\qquad
T_{ij}=U_{ij}/s_{ij},\qquad T_{ji}=T_{ij}^{-1}.
\]

Use **one common square root** for the two directions, with the branch
continuing from $+1$ for short edges. The implementation uses the principal
root of the gauge-invariant product away from its cut. It evaluates each
canonical edge only once; backward traversal uses the negative of its phase
and log magnitude.

Since $\operatorname{Log}(U_{ij}U_{ji})=-g_{\alpha\beta}\delta_\alpha\delta_\beta
+O(h^3)$, the symmetric second-order term cancels:

\[
\operatorname{Log}T_{ij}=\int_i^j a+O(h^3),\qquad
\Phi_{ijk}=i\operatorname{Log}(T_{ij}T_{jk}T_{ki}).
\]

Thus the local flux error is $O(h^3)$, or $O(h)$ after division by area for
shape-regular triangles and bounded derivatives. This estimate does not
apply uniformly at a band degeneracy or on arbitrarily thin triangles.
For Hermitian states, $U_{ji}=U_{ij}^*$ and $T_{ij}=U_{ij}/|U_{ij}|$.

Under $|R_i\rangle\mapsto c_i|R_i\rangle$ and
$\langle L_i|\mapsto c_i^{-1}\langle L_i|$, with any nonzero complex $c_i$,
$U_{ij}U_{ji}$ is unchanged and $T_{ij}\mapsto(c_j/c_i)T_{ij}$. Closed loops
are invariant. Do **not** replace this construction by independently taking
the principal root of $U_{ij}/U_{ji}$: that can introduce gauge-dependent signs.

Fukui, Hatsugai, and Suzuki use normalized Hermitian links and **inverse links
on backward traversal**, establishing a gauge-invariant lattice Chern number
and the role of admissibility [3, Eqs. (7)-(12)]. The LR square-root
normalization above is derived here; it is not an equation attributed to that
paper. The proof of continuum RR/LR equality in [1] is a separate statement.

On a closed consistently oriented triangulation, reciprocal edge factors
cancel pairwise, so the total imaginary flux vanishes and the summed principal
real flux is an integer multiple of $2\pi$, up to roundoff. This is a lattice
property. Correct recovery of the continuum integer still requires resolved
band geometry, the continuous square-root branch, and no unresolved triangle
phase winding. Integer output alone is not a convergence test.

### API, invalid inputs, and saved data

- `biorthogonal_edge_links(right, left, edges)` returns $T_{ij}$ for the supplied
  directed edges. `left` contains left **kets**; both state arrays must satisfy
  `vdot(left[i], right[i]) == 1` at each vertex.
- `triangle_flux_biorthogonal_complex(right, left, tri)` returns the full
  complex flux. Its real part uses the principal loop phase; its imaginary
  part is the signed sum of edge log magnitudes. No full loop product is
  formed, avoiding unnecessary overflow/underflow.
- `triangle_flux_biorthogonal(...)` retains the old `(real_flux, bad_count)`
  interface. Invalid LR links now raise `ValueError`; successful calls return
  `bad_count == 0`.
- `integrate_chern(verts, triangles, E, mu1, mu2, Hfun, mode="biorthogonal")` returns
  `(ChernResult, flux_data)` and records
  `link_rule="biorthogonal-reciprocal-v1"` and `total_flux_imag`. With
  the application wrapper's `save_flux=True`, the NPZ contains real `flux` and
  `flux_complex`, together with `link_rule` and oriented triangles. Right-state
  mode retains its former real flux and stores zero imaginary part.

`ChernResult` contains scalar diagnostics, without a mesh filename;
`flux_data` contains the arrays and rule metadata. The application-local
`calculate_chern` reads/writes NPZ files around this library call. The legacy
playground CLI has its own thin file/model adapter and imports the same
numerics; the application never imports that adapter.

The eigenvector routine rejects a nonfinite or numerically zero left/right
self-overlap instead of silently skipping normalization. The link routine
rejects a nonfinite/near-zero invariant product
(`abs(Uij*Uji) <= LINK_EPS**2`) and products on the negative-real square-root
cut within `LINK_BRANCH_CUT_TOL`. These checks identify an undefined or
ambiguous transport; they do not establish smoothness on a coarse edge.
Refinement/band inspection is required if an edge fails. The current method
handles a single smooth band, not a degenerate multiband subspace.

Reintegrate the saved QWZ meshes without new energy scans or refinement:

```sh
python application/Chern-number-calculation.py chern
python application/Chern-number-calculation.py diagnose
python application/Chern-number-calculation.py report
python -m pytest tests/test_chern_biorthogonal.py -q
```

Data default to `application/data/QWZ-benchmark`. The 2026-09-26 correction
archives old raw-overlap outputs under `biorthogonal-raw-overlap-baseline`;
the report checks all eight saved input mesh hashes against that archive.
Both refined LR band Chern numbers now agree with RR to roundoff: in the
application's literature convention, $(C_-,C_+)=(1,-1)$ for $m=1.2$ and
$(0,0)$ for $m=2.6$, with $\gamma=0.4$. The existing edge-length limitations
are recorded separately in the application report.

### Analytic and numerical regression checks

For $R=(1,x)^T$, $L_{\mathrm{bra}}=(1-sx,s)$, $s=1+i$, the Hamiltonian
$H=2|R\rangle\langle L|-I$ is smooth with eigenvalues $\pm1$. Its connection
is $a=s\,dx$, hence $B=0$. For vertices $(0,0),(h,0),(0,h)$ the old loop is
$1-s^2h^2=1-2ih^2$, giving apparent curvature approaching $4$ instead of
zero. Reciprocal links give zero flux.

Replacing $s$ by $(1+i)y$ gives $a=(1+i)y\,dx$ and $B=1-i$. The same
triangle now gives
$\Phi=-\tfrac{i}{2}\operatorname{Log}[1+(1+i)h^2]$; division by $h^2/2$
converges to the correct **complex** curvature. Tests also cover this example
through SVD eigenvector extraction and NPZ persistence, Hermitian reduction,
irregular QWZ tori in both phases, complex gauge changes, vertex relabeling,
orientation reversal, zero links, branch-cut products, and self-orthogonality.

## Parameters and constants

Parameter tables remain in [constants.md](constants.md). Consult its
[band-clustering settings](constants.md#pygbz2dexperimentalband_clustering)
and [mesh/Chern settings](constants.md#experimental-mesh-and-chern-modules)
for defaults. Per-call overrides, module-level assignments and multiprocessing
scope follow the [customization contract](constants.md#the-customization-model).

## Application boundary and dependency audit

The QWZ example uses the installed `pygbz2d` package, not sibling scripts.
After `python -m pip install -e .` (or installation of a current wheel),
`application/Chern-number-calculation.py` can be copied into another directory.
Its default outputs are `data/QWZ-benchmark` and
`Figures/Chern-number-calculation` beside the copied script. NumPy/SciPy are
package dependencies; BerryPy/SymPy are needed for the example's model-building
stages and Matplotlib for its reporting stage. No Python code is loaded from
`playground` or another application file.

Before this extraction, `src/pygbz2d/experimental` contained only
`band_clustering.py` and its exports. None of the mesh, refinement, or Chern
functions below was already implemented in the package.

| Category | Functions or responsibility | Location after extraction |
| --- | --- | --- |
| Already in the library | `BandPoints`, `flatten_results`, `BandClustering`, `cluster_bands` and clustering diagnostics | `experimental.band_clustering`, unchanged |
| Reusable geometry, previously only in `demo_torus_mesh.py` | `dedup_vertices`, `periodic_delaunay`, `mesh_topology`, `build_cluster_mesh` | `experimental.torus_mesh`; `build_band_mesh` also accepts an already selected point cloud |
| Reusable geometry, previously in `demo_mesh_refine.py` / `chern_on_gbz.py` | `torus_midpoint`, `torus_pair_dist`, `edge_lengths`, `triangle_areas`, `edge_len_percentiles`, `unwrap_triangle`, `orient_triangles`, `regular_torus_mesh` | `experimental.torus_mesh` |
| Reusable refinement, previously in `demo_mesh_refine.py` | `build_E_interpolator`, `chord4`, `match_result_to_target`, `dyadic_geometry_refine`, `bisect_to_boundary`, `solve_batch` and its spawn-safe worker | `experimental.mesh_refinement` |
| Reusable iteration inside `demo_pipeline.refine_to_convergence` | Repeated prediction, solving, accepting computed samples, deduplication and retriangulation | `experimental.mesh_refinement.refine_mesh`, taking arrays rather than a demo clustering object |
| Reusable Chern numerics, previously in `chern_on_gbz.py` | SVD/null-vector diagnostics, left/right eigenvectors, RR and reciprocal LR fluxes, integration | `experimental.chern`; `integrate_chern` takes arrays and returns summary plus flux arrays |
| Outer orchestration | QWZ parameters/model, energy windows and scans, band selection, exact partner-band reuse, filenames, NPZ/pickle/JSON I/O, CLI, figures and provenance | The QWZ application itself |
| Presentation only | `seam_display_filter`, topology text formatting, figure layout | Flat-plot mask implemented locally as `seam_display_mask` in the application; demo presentation remains in demos |
| Demo-only adapters | Cluster selection around `refine_to_convergence`, model-file discovery, BerryPy object adaptation, `calculate_chern` file wrapper and CLI self-test | Playground entry points, calling the package numerics |

The imported numerical implementation is unique: the affected playground
scripts now import the library routines rather than keeping independent
copies. Package modules import no playground code, model files, plotting
packages, or command-line loaders. File paths are handled by the application,
including its small `calculate_chern` NPZ wrapper. Source hashes are obtained
from loaded package modules, not fixed checkout paths; source files are not
required merely to generate a report.

## Validation and reproduction

The saved QWZ input meshes were retained. Reintegration reproduces all 16
base/refined, lower/upper, RR/LR results and their numerical diagnostics
exactly. Rebuilding all four base meshes from the saved fine-sweep subsets
also reproduces their vertex, triangle, energy and radius arrays exactly.

Regression tests cover periodic topology and seams, first-sample provenance,
live defaults, complex interpolation, insertion of solved rather than
predicted states, boundary bisection, Chern geometry/gauges, and the
application dependency boundary. A separate wheel test loads the package
from a temporary installation outside the checkout and runs the copied
example's integration, BZ checks, diagnostics and report. Two spawned worker
processes perform real GBZ root solves and match serial output.

```sh
python -m pytest tests/test_experimental.py tests/test_experimental_mesh.py tests/test_chern_biorthogonal.py tests/test_chern_example_independence.py tests/test_constants.py -q
```

## References

1. H. Shen, B. Zhen, and L. Fu, *Topological Band Theory for Non-Hermitian
   Hamiltonians*, Phys. Rev. Lett. **120**, 146402 (2018).
   [DOI](https://doi.org/10.1103/PhysRevLett.120.146402),
   [author manuscript](https://arxiv.org/abs/1706.07435).
2. Y. Ashida, Z. Gong, and M. Ueda, *Non-Hermitian Physics*, Adv. Phys. **69**,
   249-435 (2020), Sec. 5.3.
   [DOI](https://doi.org/10.1080/00018732.2021.1876991),
   [author manuscript](https://arxiv.org/abs/2006.01837).
3. T. Fukui, Y. Hatsugai, and H. Suzuki, *Chern Numbers in Discretized Brillouin
   Zone: Efficient Method of Computing (Spin) Hall Conductances*,
   J. Phys. Soc. Jpn. **74**, 1674-1677 (2005).
   [DOI](https://doi.org/10.1143/JPSJ.74.1674),
   [author manuscript](https://arxiv.org/abs/cond-mat/0503172).

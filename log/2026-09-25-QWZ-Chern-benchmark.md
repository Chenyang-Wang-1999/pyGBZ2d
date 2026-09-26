# Standard non-Hermitian QWZ benchmark

Requested workflow: 20-by-20 coarse energy scan, then 201-by-201 fine scan,
then the existing adaptive mesh refinement and Chern-number routines, with
20 workers. All numerical data are under `application/data/QWZ-benchmark`.

Model: `(sin(kx)+i*gamma)*sigma_x + (sin(ky)+i*gamma)*sigma_y
+ (m-cos(kx)-cos(ky))*sigma_z`. Cases are `gamma=0.4`, `m=1.2` and `m=2.6`.
Literature lower-band targets are C=1 and C=0, respectively (PRX 14, 021011,
Fig. 9). The existing playground code has an overall opposite Berry-phase
convention; both raw C and `literature_sign_chern=-raw_C` are saved.

## Software and method

The application is a model/driver layer. The package GBZ algorithms and
playground triangulation/refinement/Chern algorithms are used without edits.
The scan uses the NumPy Laurent backend and default solver tolerances.
BLAS threading is capped to one thread per worker to avoid oversubscription.
Source hashes for the numerical routines are saved in `source-provenance.json`.
The checkout HEAD during the run was `af595fa` (unrelated existing working-tree
changes in playground/Haldane-model-gainloss.py and paper files are preserved).

The two bands are selected by the real-energy line gap, not by fitting a
clustering radius. Every point and every stored line-subset row is retained.
The default existing refinement settings are used: edge threshold 0.2 radians,
match tolerance from `demo_mesh_refine.MATCH_TOL` (0.5), up to 10 outer rounds.
The table in the note will distinguish the base fine-sweep mesh from the
result of adaptive refinement, with no assertion that a Chern integer alone
establishes convergence.

## Completed preliminary checks

- BerryPy gives ten Laurent terms for each model. Eight deterministic complex
  arguments check the polynomial against the direct Hamiltonian determinant.
- The 20-by-20 coarse scans each completed all 400 solves without errors.
  Nontrivial: 36 in-spectrum energies, Re extent +/-3.042105263 and Im extent
  +/-0.214705627. Trivial: 16 in-spectrum energies, Re extent +/-4.294736842;
  all detected coarse spectral nodes are real.
- BZ/Hermitian checks at n=21,41,81, gamma=0 and 0.4 use the existing
  `regular_torus_mesh` and `calculate_chern` routines. Right and biorthogonal
  modes reproduce lower-band C=1 for m=1.2 and C=0 for m=2.6 to roundoff.
  They have zero bad links. These are not GBZ results.

## Coarse-window correction

The first fine window, based on two coarse-cell margins alone, was too narrow
in Im(E) for m=1.2. At its lower boundary Im(E)=-0.357842712, the solver found
nonempty subsets at Re(E)=+/-1.014631579. The calculation was stopped, with
7,000 completed entries preserved as
`QWZ-nontrivial-fine-truncated-window.pkl`. It is not used for triangulation.

This identifies a coarse-grid coverage failure: a narrow complex spectral
lobe can lie between coarse energy nodes. The corrected window uses the
finite-OBC numerical-range bound |Im(E)| <= sqrt(2)*gamma, because the only
anti-Hermitian part is the onsite gamma*(sigma_x+sigma_y). One coarse
imaginary spacing is added. Both fine grids therefore use Im(E) in
[-0.637253967,0.637253967], still with exactly 201 imaginary points.
Their real intervals retain the original two-cell margins. The script checks
for spectral points on the final outer grid before constructing the mesh.

The restarted fine sweep has again been evaluated directly at every energy;
no symmetry-generated entries and no imported Ronkin data are used.

## Completed fine sweeps

Both corrected 201-by-201 grids completed all 40,401 solves with zero failures,
zero nonfinite GBZ samples, and zero spectral hits on the outer grid.

| m | Spectral energy nodes | GBZ samples | Wall time (20 workers) | Max normalized polynomial residual |
|---|---:|---:|---:|---:|
| 1.2 | 4,302 | 25,836 | 3,651.985 s | 6.23e-14 |
| 2.6 | 1,416 | 8,440 | 1,908.422 s | 6.57e-15 |

The m=1.2 base meshes each have 12,918 vertices, 25,836 triangles, one
component, chi=0, no boundary or nonmanifold edges, and total phase area
4*pi^2. Their maximum edge is 0.666298 radians. The lower-band base Chern
numbers in the literature convention are 1.0000000000000002 (RR) and
0.9999999992982599 (LR). These are base-mesh values, before the requested
adaptive refinement.

The first nontrivial refinement round generated 2,092 targets. Prediction
energies are finite, with Re in [-3.1246,-0.9250] and Im in approximately
[-0.3734,0.3734]. The first-round positions, predicted energies and anchors
are retained in `nontrivial-first-refinement-targets.npz` for diagnostics.

Round 1 accepted all 2,092 targets: 1,798 direct, 294 by boundary fallback,
zero rejected, with 1,764 fallback probes. It took 2,144 seconds. The mesh
grew to 15,010 vertices and 30,020 triangles, and maximum edge decreased
from 0.666298 to approximately 0.3710 radians. Round 2 generated 175 targets.

For this trace-zero two-band model, f(E,beta)=f(-E,beta), so the Ronkin
function and its minimizing radii coincide at the partner energies. The
driver therefore shares the refined lower-band geometric mesh with the
upper band under E -> -E. Upper eigenvectors and Wilson phases are computed
independently. This avoids duplicating geometric refinement and is explicitly
distinct from the full energy-grid sweeps, where every node was solved.

## Final run results

The nontrivial lower-band run reached the existing ten-round cap. It has
15,274 vertices, 30,548 triangles, closed-torus topology, maximum edge
0.23192381947 rad, and eight unique edges above the requested 0.2 rad.
These remaining edges lie near the thin complex-energy tips. Later rounds
accept target matches about 0.188 rad away and increasingly merge them into
existing vertices, so additional rounds no longer reduce the maximum edge.
The target edge bound has NOT been met. The existing routines and tolerances
were left unchanged.

After refinement, lower-band C_lit is 1.0000000000000002 in right mode and
1.0058993004092458 in the existing biorthogonal triangle mode. On the shared
upper-band mesh, these values are -1.0000000000000002 and
-0.9941006995907518, respectively. All vertices have nullity one, no bad
links, and eigenvector residuals near 2e-13.

The initial process was stopped after its lower mesh and diagnostics had
been saved, avoiding a duplicate upper-band geometric refinement. The driver
was restarted to reuse that saved mesh and its exact E -> -E partner.
The trivial lower-band refinement also completed ten rounds. Round 1 accepted
1,338 of 1,352 targets (1,322 direct, 16 boundary, 14 rejected; 96 fallback
probes) in 1,204 seconds. It grew to 5,558 vertices and reduced the maximum
edge to 0.4190 radians. The final output has 5,841 vertices, 11,682 triangles,
and maximum edge 0.25348214987960016 radians, with 23 unique long edges.
It is a closed torus with no unused/nonmanifold/boundary edges, but also did
not meet the 0.2-radian target. Lower-band C_lit is 2.2087185287941088e-16 (RR)
and -1.874131958518914e-5 (LR). The shared upper-band values are zero (RR)
and -1.8741319579226705e-5 (LR).

All requested computational stages finished on 2026-09-26 local time.
The default right-state benchmark reproduces the reference integers. The
global edge-length target and experimental LR discretization do not pass
their respective checks. These limitations are retained explicitly rather
than reported as a fully converged, all-modes benchmark.

## Biorthogonal discretization diagnosis

The discrepancy is reproducible without any further GBZ solve. Reverse every
triangle on the same saved nontrivial lower-band mesh:

| Existing mode | Raw C, positive orientation | Raw C, negative orientation |
|---|---:|---:|
| Right | -1.0000000000000002 | 0.9999999999999999 |
| Biorthogonal | -1.0058993004092458 | 0.9941006995907546 |

The existing LR rule uses the three independent overlaps around each triangle,
U_ij=<L_i|R_j>. Generally U_ji is not the conjugate/inverse of U_ij. Therefore
the two faces sharing an edge contribute a nonzero phase from U_ij*U_ji.
For this mesh, direct evaluation gives

`-sum_unique_edges Arg(U_ij*U_ji)/(2*pi) = -0.0058993004092451826`.

This accounts for the entire fractional offset of the raw LR result and for
the nonzero even component under orientation reversal. It isolates a defect
of the experimental finite-mesh LR triangle rule, rather than a failure of
the GBZ eigenvalue residual or of the continuous biorthogonal Chern definition.
The original symmetric base mesh happens to cancel this contribution nearly
exactly; that accidental cancellation is lost after irregular refinement.
No correction or rounding has been applied in the benchmark tables.

The reproducible `diagnose` command saves the orientation checks, shared-edge
phase contribution and remaining long-edge endpoints. Repair of the existing
LR integrator is a separate code modification; a consistent discrete link must
respect edge reversal and be checked on non-Hermitian irregular meshes.

For the trivial refined mesh, the LR shared-edge bias is
`+1.874131958389269e-5`. Both raw orientation choices give approximately this
same positive value, while RR is zero to roundoff. This independently
reproduces the same LR discretization defect in the C=0 case.

## Delivery validation

- All stages (`coarse`, `fine`, `mesh`, `bz-check`, `diagnose`, `report`) were
  executed on the saved benchmark data.
- Both refined tori pass the closed-manifold checks; the global 0.2-radian
  edge target remains unmet in both cases and is reported as such.
- The recorded SHA-256 hashes of all package/playground numerical routines
  remain unchanged after the run.
- The application script passes Python bytecode compilation. Changed tracked
  text files pass `git diff --check`.
- The final seven-page Typst PDF was compiled and every page visually checked.
  Fine-sweep panels explicitly indicate their vertical zoom; Berry-flux maps
  use a common, zero-centered signed color scale. Original numerical values,
  including the LR defects and remaining long edges, are retained.

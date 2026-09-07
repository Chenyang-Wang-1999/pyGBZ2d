# Sources and scope

## Theory and derivation

- Hong-Yi Wang, Fei Song, and Zhong Wang, *Amoeba Formulation of Non-Bloch Band Theory in Arbitrary Dimensions*, Physical Review X 14, 021011 (2024). [DOI](https://doi.org/10.1103/PhysRevX.14.021011), [full text](https://arxiv.org/html/2212.11743v3). Ronkin minimum and central-hole interpretation: PDF page 3. Pages 7–8 explicitly derive the root-modulus and root-count representation by polynomial factorization and Jensen's identity.
- Chenyang Wang et al., *Universal theory for geometry-dependent non-Hermitian bands*, [arXiv:2506.22743v2](https://arxiv.org/abs/2506.22743v2), [full text](https://arxiv.org/html/2506.22743v2). Strip physical picture and its distinction from the amoeba construction: pages 2, 4–5. Transverse middle-root constraint: page 9. Equivalent loop-winding formulation (Supplement S2): page 10. The explicit version is retained because the unversioned manuscript title may differ.

The project contribution presented here is a numerical bridge and its design, not a claim to have originated the two physical GBZ formulations or Jensen's identity. A common root representation does not identify the two GBZs.

## Assumptions behind the bridge

- Use beta_j = exp(mu_j + i theta_j). With this convention, a wave proportional to beta_j^n has envelope exp(mu_j n).
- f = det(EI-h) is a finite Laurent polynomial. M is its transverse pole order; the generic number of transverse roots is M+N, counted with multiplicity. The displayed factorization assumes a nonzero leading coefficient and finite nonzero roots. Degree loss and multiple roots require limiting or separate treatment.
- Jensen's identity retains ln|c_N| and the pole term -M mu2. Neither may be discarded from the full Ronkin expression. The leading coefficient is constant with respect to mu2 at fixed beta1, but generally varies with mu1.
- The root-count derivative is stated at differentiable points. Persistent root-level equality requires one-sided counts or the convex subgradient condition. Exchanging axes gives the mathematical a1 identity; this does not describe a separate axis-swapped computation in the implementation.
- The sorted middle-root gap is nonnegative. Locating it numerically cannot rely only on a sign change of that sorted gap. Tracked-root differences and event geometry resolve candidate ties.
- Loop winding uses the total derivative along the varying middle-gap reference. Its theta2 average can be noninteger. Singular loops require one-sided treatment. Zero-free plateaus are classified separately from supported GBZ subsets.
- The separable Hatano–Nelson benchmark on page 16 follows analytically from the displayed Hamiltonian, real positive J and real gamma, with axis-aligned OBC. It supplies no evidence for generic geometry-dependent accuracy or runtime.

## Current-code design map

The code takes precedence over older manuscript and pseudocode descriptions. This map keeps implementation references available without making them the talk's main narrative.

| PDF pages | Principal repository sources |
| --- | --- |
| 6 | `src/pygbz2d/core.py`, `backend.py` |
| 7–8 | `src/pygbz2d/amoeba/ronkin_winding.py`, `zm_extract.py`, `bisect.py` |
| 9–11 | `src/pygbz2d/sgbz/mu2mid.py`, `pairwise.py`, `winding.py`, `sgbz_solver.py`, `plateau.py` |
| 12–14 | `src/pygbz2d/continuation/zero_manager.py`, `arclength.py`, `multiple_roots.py` |
| 15 | `src/pygbz2d/core.py`, `amoeba/amoeba.py`, `sgbz/continuum_lines.py` |
| 16 | `conftest.py::build_HN2D_polynomial` in the explicitly stated separable limit |
| 17 | `DATA_REQUIREMENTS.md`, `collect_slide_data.py`, `plot_slide_data.py`: proposed validation, no empirical claim |
| 18 | Experimental band reconstruction; `experimental/band_clustering.py`, `playground/demo_pipeline.py` |

Unresolved close events, exact degeneracies and finite sampling remain numerical limitations. No speedup, completed benchmark, or universal numerical certification is claimed. Warm starts across radius probes remain proposed work and are not presented as implemented.

## Artwork and typography

The editable SVG schematics are conceptual geometry and design diagrams, with matching vector PDFs. None contains fabricated sampled model data. Three schematics appear in the main deck; the nested amoeba-search diagram remains available as a supplementary asset. The pinned Touying 0.7.4 University theme retains the user's academic style. Equations and text remain editable in Typst.

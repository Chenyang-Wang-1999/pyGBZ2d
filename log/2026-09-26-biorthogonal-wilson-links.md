# 2026-09-26: reciprocal biorthogonal Wilson links

The user approved the theoretical correction and requested source/doc updates
with references. This supersedes the unresolved LR-link diagnosis in
`2026-09-25-QWZ-Chern-benchmark.md`; the earlier numerical record is retained.

## Change and derivation

Raw links `Uij = <Li|Rj>` carry a second-order symmetric quantum-geometric
contribution. Its imaginary part contaminates the phase at the same order as
a shrinking triangle's physical flux. Raw `Uij*Uji` also contributes a phase
on internal edges, explaining the orientation-even offset of the old rule.

`playground/chern_on_gbz.py` now uses `Tij = Uij/sqrt(Uij*Uji)`, computing one
common square root per undirected edge and imposing inverse backward
transport. The complex flux is `i Log(Tij*Tjk*Tki)`. Signed edge phase and
log-magnitude sums avoid full loop products. The old real-flux helper API is
retained; an additional helper returns complex flux. Saved NPZs keep `flux`
and add `flux_complex` and the rule tag `biorthogonal-reciprocal-v1`. Result
JSON adds `link_rule` and `total_flux_imag`.

Zero/ill-conditioned self-overlaps are explicitly rejected, as are undefined
edge products and products on the square-root branch cut. A diagnostic fix
also avoids NaN percentile interpolation when every singular-value gap ratio
is infinite (encountered in the exact analytic regression model).

The derivation and validity conditions are in `doc/chern.md` and the updated
application Typst/PDF. Citations distinguish the RR/LR equality theorem of
Shen-Zhen-Fu, the definitions reviewed by Ashida-Gong-Ueda, and inverse
backward links/admissibility in Fukui-Hatsugai-Suzuki. The particular LR
normalization is derived here rather than attributed as an original FHS
formula. No novelty claim is made for it.

## Preserved inputs and recalculated results

Old Chern/diagnostic JSON, biorthogonal flux archives, BZ checks, integrator
source, summary, and source hashes were copied to
`application/data/QWZ-benchmark/biorthogonal-raw-overlap-baseline` before
recalculation. Its input-hash manifest covers all eight saved base/refined
band meshes. The new `chern` application stage reintegrates those files
without rewriting them. No GBZ scan or adaptive refinement was rerun.

Executed `chern`, `bz-check`, `diagnose`, and `report`. All eight mesh hashes
match the archive. In the application's literature sign convention:

| mass | lower RR / LR | upper RR / LR |
| --- | --- | --- |
| 1.2 | +1 / +1 | -1 / -1 |
| 2.6 | 0 / 0 | 0 / 0 |

Agreement is at floating-point accuracy. The legacy lower-band raw-sign
shared-edge biases are still recorded as historical diagnostics:
`-0.0058993004092451826` (nontrivial) and
`+1.874131958389269e-5` (trivial). They are not residual errors in the
corrected rule.

On refined lower-band meshes, the largest local complex-gauge discrepancy
is `5.8364e-15`, local orientation reversal error is zero, and independent
forward/reverse edge products differ from one by at most `8.8819e-16`.
The imaginary integrated Chern values are below `3e-18`. Invariant edge
products have minimum real parts `0.9869448554` and `0.9709417958`.

The existing long-edge counts (8 and 23), ten-round refinement cap, and
maximum lengths remain unchanged. Integer lattice Chern output is not used
to certify geometric or local-curvature convergence.

## Regression verification

`python -m pytest tests/test_chern_biorthogonal.py -q`: **10 passed**.
Coverage includes a flat connection whose old apparent curvature approaches
4, an analytic complex curvature `1-i`, Hermitian reduction, both phases of
non-Hermitian QWZ on irregular tori, arbitrary complex gauges, vertex
relabeling, orientation, zero/branch-cut links, exceptional self-overlap,
SVD eigenvectors, and persistence of complex flux metadata.

The QWZ BZ checks use 21, 41, and 81 points per direction, at gamma 0 and
0.4. The application report preserves all unrounded values and current
source hashes, and the updated Typst source compiles to PDF.

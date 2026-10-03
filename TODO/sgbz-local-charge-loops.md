# TODO: recover unknown SGBZ charges with local loops

Status: deferred (2026-10-03). The chosen immediate fix is to omit sufficiently
narrow theta2 regions using a tolerance tied to the winding calculation.

## Motivation

`diagnostics/Haldane-y-SGBZ.py` reproduces a failure at `E=-1.4675` for the
y supercell. At the first plateau probe, `mu1=-1e-7`, MR events with unknown
charges delimit a theta2 region of width approximately `2.56635e-7` across
the periodic seam. Without a charge to propagate, the winding algorithm
must choose a seed inside this narrow region.

The alternative proposed by the user is to recover an unknown charge from
a loop enclosing its zero. This would allow winding propagation through
more events and reduce the need for separate seed integrations.

## Proposed calculation

For fixed `(E, mu1)`, use the same function as the major-axis winding:

```
G(theta1, theta2) = f(E, exp(mu1 + i*theta1),
                      exp(mu2_mid(theta1) + i*theta2))
```

Compute the change in `arg(G)` on a closed contour around an isolated zero
in the `(theta1, theta2)` torus. Calibrate contour orientation against a
known ordinary event: the stored charge is the jump in the major-axis loop
winding as theta2 increases, so the sign convention must be checked rather
than assumed from a counterclockwise contour.

## Open issues

- Radius selection must account for both angular coordinates and periodic
  seams. Proximity of theta2 projections alone does not imply that zeros
  are close on the torus.
- A small radius needs a discretization-error check: refine contour samples,
  resolve phase increments, and verify the recovered integer is stable.
  An apparently integer answer from one coarse polygon is insufficient.
- A large radius may enclose other zeros. Bound it using neighbouring zeros
  and verify that the contour itself does not intersect a zero. Account for
  uncertainty in the zero locations.
- Coincident zeros may permit only a combined charge to be measured. Define
  how a combined charge enters propagation without assigning arbitrary
  individual charges to snapped root columns.
- The contour uses the piecewise mu2_mid path, including MR boundaries and
  derivative jumps. Distinguish error in that path and in event locations
  from error in sampling the enclosing contour.
- Compare results across a usable range of radii and contour resolutions.
  If isolation and charge convergence cannot be established, retain the
  unknown charge instead of silently trusting the estimate.

## Acceptance cases

- Ordinary crossings with known positive, negative, and zero charge.
- Isolated MRs, tangencies, and nearly coincident zeros.
- Zeros across either periodic seam and zeros with equal theta2 but distant
  theta1 coordinates.
- A contour enclosing multiple zeros, checked against their total charge.
- The Haldane-y reproducer and agreement with independent loop windings in
  well-resolved neighbouring theta2 regions.

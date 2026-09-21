#import "@preview/zhaji:0.1.0": * 

#let nt = note(
  title: "Benchmark: 2D Hatano-Nelson model",
  author: "Chenyang Wang",
  font-head: "Arial",
	font-text: "Times New Roman",
	first-line-indent: 0em,
	lang: "en"
)

#show: nt.make
#show raw.where(block: true): it => block(it, fill: rgb("edf6fd"), inset: 5pt, breakable: false)
#show raw.where(block: false): it => highlight(it, fill: rgb("d0eafb"), top-edge: 1.2em, bottom-edge: -.5em, extent: 0.2em)

#let mathbf(x) = $bold(upright(#x))$
#let rme = $upright(e)$
#let rmi = $upright(i)$


= Benchmark: 2D Hatano-Nelson model

This example benchmarks `pyGBZ2d` using the two-dimensional Hatano-Nelson (HN) model. Closed-form expressions for its amoebic generalized Brillouin zone (GBZ) and strip GBZs provide a reference for checking the numerical results @wang2025generaltheorygeometrydependentnonhermitian (see #link("https://arxiv.org/abs/2506.22743v3")[arXiv:2506.22743v3]).

The accompanying #link("./benchmark-2D-Hatano-Nelson.py")[Python script] demonstrates how to:
- Construct a characteristic polynomial from coefficients and exponents.
- Compute amoebic and strip-GBZ subsets at a reference energy $E_("ref")$.
- Read and visualize `PointSubset` and `LineSubset` data.
- Compare every returned sample with the closed-form GBZ and check its energy residual.
- Check whether an empty or nonempty result agrees with the analytic spectrum.

Run the script from the repository root with `python application/benchmark-2D-Hatano-Nelson.py`. It uses the local `src` directory by default and requires NumPy, SciPy, and Matplotlib for plotting. Add `--compare-only` to run the numerical checks without importing Matplotlib or opening figures. Set `RUN_IN_SRC = False` to use an installed copy of `pyGBZ2d` instead.

== Introduction

The Hamiltonian of the 2D HN model is
$ H = sum_(i , j) (J_(x 1) c_(i + 1 , j)^(dagger) c_(i , j) + J_(x 2) c_(i - 1 , j)^(dagger) c_(i , j) + J_(y 1) c_(i , j + 1)^(dagger) c_(i , j) + J_(y 2) c_(i , j - 1)^(dagger) c_(i , j)) , $
where $c_(i , j)$ annihilates a particle at site $(i,j)$ and the hopping amplitudes are complex. Throughout this benchmark, all four hopping amplitudes are nonzero. Each pair can be written as
$ J_(alpha 1) = rme^(gamma_(alpha) + rmi delta_(alpha)) J_(alpha) , quad J_(alpha 2) = rme^(- gamma_(alpha) + rmi delta_(alpha)) J_(alpha)^(*) , quad alpha = x , y $
where $gamma_(alpha), delta_(alpha) in RR$ describe the hopping asymmetry and common phase, respectively. The pair $J_(alpha), J_(alpha)^*$ forms a Hermitian hopping pair before these factors are applied. A branch convention for this decomposition is
$ gamma_(alpha) = & frac(1, 2) ln (frac(abs(J_(alpha 1)), abs(J_(alpha 2)))) \
delta_(alpha) = & frac(1, 2) limits("Arg") (J_(alpha 1) J_(alpha 2)) \
J_(alpha) = & J_(alpha 1) \/ rme^(gamma_(alpha) + rmi delta_(alpha)) $

Here, $limits("Arg")$ denotes the principal argument, as implemented by `numpy.angle`. The decomposition has a sign ambiguity: replacing $(delta_(alpha), J_(alpha))$ by $(delta_(alpha)+pi, -J_(alpha))$ leaves both hopping amplitudes unchanged. The convention above fixes a representative. The script implements it in `factorize_coefficients` and checks reconstruction in `check_factorization`.

The $x$-strip, $y$-strip, and amoebic GBZs coincide when expressed in the same Cartesian coordinates:
$ beta_(x) & = exp (gamma_(x) + rmi theta_(x)) , \
beta_(y) & = exp (gamma_(y) + rmi theta_(y)) , $
where $theta_x, theta_y in [0, 2 pi)$. The $[11]$-strip GBZ generally differs. In the basis $(mathbf(a)_([11]), mathbf(a)_y)$, with $mathbf(a)_([11]) = mathbf(a)_x + mathbf(a)_y$, it is
$ tilde(beta)_([11]) & = rme^(gamma_(x) + gamma_(y) + rmi theta_([11])) , \
tilde(beta)_(y) & = rme^(gamma_(y) + rmi theta_(y)) sqrt(lr(abs(frac(J_(x)^(*) rme^(rmi Delta_(x y) + rmi theta_([11])) + J_(y), J_(x) rme^(rmi Delta_(x y) - rmi theta_([11])) + J_(y)^(*))))) , $
where $Delta_(x y) := delta_(x) - delta_(y)$. These expressions follow Eqs. (S3.16) and (S3.30) of Ref. @wang2025generaltheorygeometrydependentnonhermitian. The second expression applies where its numerator and denominator are nonzero. When $sin(Delta_(x y)) eq.not 0$, the transverse radius depends on $theta_([11])$, and the model exhibits geometry-dependent bands. When $sin(Delta_(x y)) = 0$, the radius reduces to $exp(gamma_y)$ at nonsingular phases.

The corresponding energy spectra are illustrated in #ref(<fig:HN-model>).

#figure(
	image("Figures/HN-example-20250524.pdf"),
	caption: [ 
  Strip geometries and GBZ spectra of the 2D HN model. *a* The $x$-, $y$-, and $[11]$-strip geometries. *b* The common spectrum of the $x$-strip, $y$-strip, and amoebic GBZs. *c* The $[11]$-strip GBZ spectrum.
]
)<fig:HN-model>

== Modeling with characteristic polynomials
`pyGBZ2d` represents a tight-binding model through its characteristic Laurent polynomial $det(E I - h(beta_1, beta_2))$. Write its expansion as
$ f (E , beta_(1) , beta_(2)) = sum_(j = 1)^(N_("terms")) a_(j) E^(l_(j)) beta_(1)^(m_(j)) beta_(2)^(n_(j)) , quad l_(j) in NN , m_(j) , n_(j) in ZZ $
where $l_j$ is a nonnegative integer. The coefficient array has shape `(N_terms,)`, and the integer degree array has shape `(N_terms, 3)`. Their entries correspond term by term. Schematically,
```python
from pygbz2d import core
import numpy as np

coeffs = np.array([a1, a2, ..., aN], dtype=complex)
degs = np.array([
	[l1, m1, n1],
	[l2, m2, n2],
	...,
	[lN, mN, nN],
], dtype=int)
poly = core.CharPoly(coeffs, degs)
```

Changing the lattice basis changes the coordinates used to describe a GBZ. Under an invertible integer change of lattice basis, the amoebic GBZ describes the same physical set after transforming the coordinates back. Strip GBZs additionally depend on the choice of major axis, so the two input Bloch factors play different roles.

The `pygbz2d.sgbz` module assumes that the polynomial is #text(c-emph)[expressed in the basis $(mathbf(a)_1, mathbf(a)_2)$, with $mathbf(a)_1$ as the major axis]. A change of major axis can therefore change the physical GBZ, even after both results are expressed in Cartesian coordinates.

For the Cartesian basis $(mathbf(a)_x, mathbf(a)_y)$, the Bloch Hamiltonian is

$ h (beta_(x) , beta_(y)) = J_(x 1) beta_(x)^(- 1) + J_(x 2) beta_(x) + J_(y 1) beta_(y)^(- 1) + J_(y 2) beta_(y) , $
and the characteristic polynomial is
$ f (E , beta_(x) , beta_(y)) = E - J_(x 1) beta_(x)^(- 1) - J_(x 2) beta_(x) - J_(y 1) beta_(y)^(- 1) - J_(y 2) beta_(y) . $

For a $y$-strip calculation, interchange the last two columns of the degree array so that $(beta_1,beta_2)=(beta_y,beta_x)$.

For $(mathbf(a)_([11]), mathbf(a)_y)$, the coordinate transformation is
$ beta_([11]) = beta_x beta_y, quad beta_x = beta_([11]) / beta_y. $
Substitution gives
$ f (E , beta_([11]) , beta_(y)) = E - J_(x 1) beta_([11])^(- 1) beta_(y) - J_(x 2) beta_([11]) beta_(y)^(- 1) - J_(y 1) beta_(y)^(- 1) - J_(y 2) beta_(y) . $

The function `get_HN_charpoly` implements all three choices using `which="x-y"`, `"y-x"`, or `"11-y"`. It returns `(coeffs, degs)`, which can be passed directly to either solver; constructing a `CharPoly` explicitly is optional for these entry points.

For larger models, entering coefficient and degree arrays manually can become cumbersome. The `TightBinding` module of #link("https://github.com/Chenyang-Wang-1999/BerryPy")[BerryPy] offers an alternative interface based on hopping terms. The present benchmark uses explicit arrays to make the coordinate conventions transparent and requires no additional modeling package.

== Computing GBZ subsets at a reference energy
Both `sgbz` and `amoeba` provide `collect_GBZ_subsets` to compute the intersection of the GBZ with $f(E_("ref"),beta_1,beta_2)=0$. This is the non-Hermitian counterpart of an equal-energy contour, although the intersection may consist of isolated points. The functions share the basic inputs `coeffs`, `degs`, and `E_ref`, and both return a `GBZResult`:
```python
from pygbz2d import amoeba, sgbz

coeffs, degs = get_HN_charpoly(
    1 + 1j, 1.5 + 1.2j, -1 + 1j, -1.2 - 0.5j, "x-y",
)
res_amoeba = amoeba.collect_GBZ_subsets(coeffs, degs, 1 + 1j)
res_xstrip = sgbz.collect_GBZ_subsets(coeffs, degs, 1 + 1j)
```

Their optional parameters differ. For example, `sgbz` accepts `mu1_guess`, whereas `amoeba` accepts `mu1_low`, `mu1_high`, `mu2_low`, and `mu2_high`. Both accept `perc`, an optional progress fraction printed as a percentage, and `debug_mode=True`, which propagates solver exceptions for diagnosis. All numerical options are explicitly named; there is no catch-all keyword dictionary. See #link("../doc/constants.md")[`doc/constants.md`] for the live default settings.

=== Reading the result
`GBZResult`, defined in `src/pygbz2d/core.py`, stores the reference energy, solver status, and connected subsets:

#table(
  columns: (auto, 1fr),
  inset: 6pt,
  table.header([*Attribute*], [*Meaning*]),
  [`E_ref`], [Reference energy supplied to the solver.],
  [`success`, `error`], [Whether the computation completed without error, and an error message if it did not.],
  [`subsets`], [A list of `PointSubset` and/or `LineSubset` objects.],
  [`index`], [The pair $(n_("0D"), n_("1D"))$ counting point and line subsets.],
  [`is_empty`], [Whether `subsets` is empty, regardless of solver status.],
  [`is_gbz`], [Whether `success` is true and `index` differs from `(0, 0)`.],
)

Always check `success` before interpreting an empty result. A failed calculation does not establish that the energy lies outside the spectrum:
```python
res = res_xstrip
if not res.success:
    raise RuntimeError(res.error)
if res.is_empty:
    print("No GBZ subsets found at this reference energy.")
else:
    print(res.index)  # (number of point subsets, number of line subsets)
```

The subset data fields are summarized below; validation methods and convenience properties are omitted:
```python
@dataclass(frozen=True)
class PointSubset:
    """An isolated point in the solver's coordinate basis."""
    E: complex
    beta1: complex
    beta2: complex

@dataclass
class LineSubset:
    """One sampled curve with constant mu1 = ln|beta1|."""
    E: complex
    mu1: float
    theta1_arr: np.ndarray  # shape (N,), ordered phase samples
    beta2_arr: np.ndarray   # shape (N,), eagerly stored Bloch factors
```

For a line subset, reconstruct the first Bloch factor using
```python
beta1_arr = np.exp(subset.mu1 + 1j * subset.theta1_arr)
beta2_arr = subset.beta2_arr
```
The number of stored samples is distinct from the number of line subsets. Each `LineSubset` represents one $beta_2$ curve; a degeneracy can give multiple curve objects. `subset_samples` provides a common array interface for points and lines, and `calculate_subset_and_show` demonstrates solving, comparison, data retrieval, and visualization. Phases are plotted modulo $2 pi$ using dots to avoid artificial connections across the angular seam. The dashed radius curve describes the full analytic GBZ constraint, rather than only the points at $E_("ref")$.

=== Point and line examples
For the Cartesian GBZs, substituting the analytic radii into the Hamiltonian gives
$ E(theta_x,theta_y) = 2 rme^(rmi delta_x) limits("Re")(J_x rme^(-rmi theta_x)) + 2 rme^(rmi delta_y) limits("Re")(J_y rme^(-rmi theta_y)). $
Their spectrum has nonzero area when $sin(Delta_(x y)) eq.not 0$. At a regular energy in a two-dimensional spectral region, fixing a complex energy generally imposes two real constraints and leaves isolated GBZ points. When the spectrum collapses to a line, regular energies instead generally have continuous preimages. Boundary energies and degeneracies require separate consideration.

This area statement must be qualified for the $[11]$ strip: when $abs(J_x)=abs(J_y)$, its spectrum can have zero area even if $sin(Delta_(x y)) eq.not 0$, while the GBZ remains geometry dependent. See Sec. S6 of Ref. @wang2025generaltheorygeometrydependentnonhermitian.

The script uses two fixed parameter sets, each evaluated for all four GBZ choices:
- `demo_point_subsets`: $J_(x 1)=1+rmi$, $J_(x 2)=1.5+1.2rmi$, $J_(y 1)=-1+rmi$, $J_(y 2)=-1.2-0.5rmi$, and $E_("ref")=1+rmi$. The hopping parameters match the example in Ref. @wang2025generaltheorygeometrydependentnonhermitian.
- `demo_line_subsets`: $J_(x 1)=1$, $J_(x 2)=1.5$, $J_(y 1)=-1$, $J_(y 2)=-1.2$, and $E_("ref")=1$. Here $delta_x=delta_y=0$ in the chosen convention, and the spectra are real.

== Comparison to closed-form solutions
The comparison first checks whether the reference energy belongs to the analytic GBZ spectrum, including when the numerical result is empty. For a nonempty result, it then evaluates the analytic radius constraints at every returned phase and independently checks $E_("ref")=h(beta_x,beta_y)$. Both sample checks are needed: satisfying the characteristic equation alone does not establish GBZ membership, and satisfying the radius constraints alone does not fix the energy.

=== Checking spectral membership, including empty results
Define the complex spectral half-widths
$ u=2 abs(J_x) rme^(rmi delta_x), quad v=2 abs(J_y) rme^(rmi delta_y). $
For the amoebic, $x$-strip, and $y$-strip GBZs, the dispersion relation gives
$ sigma_("Cartesian") = {u s+v t : s,t in [-1,1]}. $
Thus the spectrum is the closed parallelogram with vertices $u+v$, $-u+v$, $-u-v$, and $u-v$. When $u$ and $v$ are collinear, this set reduces to a line segment. A point is outside if it lies beyond the segment endpoints or off the segment's supporting line.

For the $[11]$ strip, write $h=A beta_2+B/beta_2$, with $A=J_(x 1)/beta_1+J_(y 2)$ and $B=J_(x 2) beta_1+J_(y 1)$. On its analytic GBZ, the equal-modulus condition gives $E=2 sqrt(A B) s$ for $s in [-1,1]$. As the phase of $beta_1$ varies, $4 A B$ traverses the line segment joining $(u+v)^2$ and $(u-v)^2$. Squaring the energy therefore gives the following equivalent membership test:
$ E in sigma_([11]) quad <==> quad E^2 in limits("conv"){0,(u+v)^2,(u-v)^2}, $
where $limits("conv")$ denotes the convex hull. This is a triangle in the squared-energy plane, possibly collapsed to a segment. Both signs of $E$ occur, so squaring introduces no extra spectral branch. The test follows directly from the strip dispersion in Eq. (S3.31) of Ref. @wang2025generaltheorygeometrydependentnonhermitian.

`closed_form_spectrum_membership` evaluates these sets independently of the numerical solver. It computes the distance to the closed polygon using edge projections and an interior test, including collapsed polygons; no angular sampling mesh is needed. Let $L=abs(u)+abs(v)$. The distance is normalized by $L$ in the Cartesian energy plane and by $L^2$ in the $[11]$ squared-energy plane. A normalized distance no larger than `spectrum_atol = 1e-8` is classified as inside. Points within this tolerance of a boundary are consequently treated as spectral points; the squared-energy tolerance has a different interpretation near $E=0$ from a tolerance applied directly to $E$.

`compare_to_closed_form` compares this analytic classification with the numerical result:
#table(
  columns: (auto, auto, 1fr),
  inset: 6pt,
  table.header([*Analytic classification*], [*Numerical result*], [*Comparison*]),
  [Inside], [Nonempty], [Proceed to the radius and energy checks.],
  [Outside], [Empty], [Pass the spectral-membership check; there are no samples to check.],
  [Inside], [Empty], [Raise a spectrum-membership mismatch.],
  [Outside], [Nonempty], [Raise a spectrum-membership mismatch.],
)
A failed solve is always an error, regardless of the analytic classification. For a verified empty result, `n_samples` is zero and the three sample-error fields are `None` (printed as `N/A`), rather than zero. The report includes `analytic_is_gbz`, `numerical_is_gbz`, `spectrum_distance`, and `spectrum_energy_power`; the last field is 1 for a Cartesian spectrum and 2 for the $[11]$ spectrum.

=== Coordinate conventions and analytic radii
Let $mu_j=ln abs(beta_j)$ and $theta_1=limits("Arg")(beta_1)$. `closed_form_log_radii` returns the following analytic values in the solver's own basis:

#table(
  columns: (auto, auto, auto, auto),
  inset: 6pt,
  table.header([*GBZ choice*], [*$(beta_1,beta_2)$*], [*$mu_1^"exact"$*], [*$mu_2^"exact"$*]),
  [amoeba / x-strip], [$(beta_x,beta_y)$], [$gamma_x$], [$gamma_y$],
  [y-strip], [$(beta_y,beta_x)$], [$gamma_y$], [$gamma_x$],
  [11-strip], [$(beta_x beta_y,beta_y)$], [$gamma_x+gamma_y$], [$gamma_y + (1/2) ln abs(R(theta_1))$],
)

Here,
$ R(theta_1) = frac(J_x^* rme^(rmi(Delta_(x y)+theta_1))+J_y, J_x rme^(rmi(Delta_(x y)-theta_1))+J_y^*). $
The $[11]$ formula can also be checked directly: in this basis the Hamiltonian is $h=A beta_2+B/beta_2$, where
$ A=J_(x 1)/beta_1+J_(y 2), quad B=J_(x 2) beta_1+J_(y 1). $
The two roots of $A beta_2^2-E beta_2+B=0$ have product $B/A$, so their equal-modulus condition gives $abs(beta_2)=sqrt(abs(B/A))$. At $abs(beta_1)=exp(gamma_x+gamma_y)$, this reproduces the radius above. The benchmark requires finite, nonzero radii and reports singular expressions rather than clipping or replacing them. The two demonstration parameter sets have $abs(J_x) eq.not abs(J_y)$, so neither the numerator nor the denominator of $R$ vanishes for real $theta_1$.

=== Error measures
For the $N$ returned samples, define
$ epsilon_(mu j) = max_(1 <= n <= N) abs(ln abs(beta_(j,n)) - mu_j^"exact"(theta_(1,n))), quad j=1,2. $
For a `PointSubset`, there is one sample. For a `LineSubset`, all stored samples are used, with $beta_(1,n)=exp(mu_1+rmi theta_(1,n))$.

To check the energy, first convert the samples to Cartesian coordinates: interchange the two components for a $y$ strip, and use $beta_x=beta_1/beta_2$, $beta_y=beta_2$ for a $[11]$ strip. Define the four hopping contributions as
$ (t_(1,n),t_(2,n),t_(3,n),t_(4,n)) = (J_(x 1)/beta_(x,n), J_(x 2) beta_(x,n), J_(y 1)/beta_(y,n), J_(y 2) beta_(y,n)). $
The dimensionless energy residual is
$ epsilon_E = max_(1 <= n <= N) frac(abs(E_("ref")-sum_(q=1)^4 t_(q,n)), abs(E_("ref"))+sum_(q=1)^4 abs(t_(q,n))). $
This normalization uses the magnitudes of the individual hopping terms, so cancellation in their sum does not create an artificially small denominator. It also remains well defined at $E_("ref")=0$ for the nonzero finite hoppings and Bloch factors considered here.

For a nonempty result with matching spectral membership, `compare_to_closed_form` returns the subset counts, sample count, and these three maximum errors. It rejects invalid samples and accepts the sample comparison when $epsilon_(mu 1),epsilon_(mu 2) <= 10^(-5)$ and $epsilon_E <= 10^(-7)$. These are benchmark acceptance thresholds, separate from both `spectrum_atol` and the solver settings; they do not modify solver tolerances. Exceeding either threshold raises an exception with the measured errors. The point and line demos also require nonempty results of their expected subset type.

=== Running the comparison
From the repository root, run
```sh
python application/benchmark-2D-Hatano-Nelson.py --compare-only
```
This executes the factorization reconstruction check, the eight point/line comparisons, and sixteen additional spectral-membership comparisons. To inspect one result from within the accompanying script or an interactive session with its functions loaded, use
```python
hoppings = (1 + 1j, 1.5 + 1.2j, -1 + 1j, -1.2 - 0.5j)
res = calculate_subsets(1 + 1j, *hoppings, "11-strip")
report = compare_to_closed_form(res, *hoppings, "11-strip")
print_comparison(report)
```

`demo_spectrum_membership` covers the following additional reference energies for all four GBZ choices. “Complex” and “real” refer to the hopping parameters in `demo_point_subsets` and `demo_line_subsets`, respectively:
#table(
  columns: (auto, auto, 1fr, 1fr),
  inset: 6pt,
  table.header([*Hoppings*], [*$E_("ref")$*], [*amoeba / x / y*], [*11-strip*]),
  [Complex], [$6+6rmi$], [Outside], [Outside],
  [Real], [$6$], [Outside: beyond endpoints], [Outside: beyond endpoints],
  [Real], [$rmi$], [Outside: off the real axis], [Outside: off the real axis],
  [Complex], [$2+2rmi$], [Inside], [Outside],
)
All sixteen additional comparisons passed in the current checkout: thirteen returned verified empty results and three returned nonempty point subsets. The last case checks geometry-dependent membership explicitly. For example, replacing `1 + 1j` by `2 + 2j` in the code above gives a verified empty $[11]$ result, with `index=(0, 0)`, `n_samples=0`, and `analytic_is_gbz=False`.

The optional `--random-benchmark --compare-only` mode uses random hoppings and samples an energy from the Cartesian GBZ. That guarantees Cartesian spectral membership, but does not guarantee $[11]$-strip membership. Its results are therefore judged by the same analytic membership test, rather than by requiring every solver result to be nonempty.

The following results were obtained on 2026-09-21 with the current checkout, the NumPy polynomial backend (NumPy 1.26.4, SciPy 1.11.4), and the default solver settings. $N$ counts all stored samples, including any shared endpoints. Adaptive sampling and future solver changes can alter $N$ and the last digits of the errors.

#figure(
  table(
    columns: (auto, auto, auto, auto, auto, auto),
    inset: 5pt,
    table.header([*Case / GBZ*], [*`index`*], [*$N$*], [*$epsilon_(mu 1)$*], [*$epsilon_(mu 2)$*], [*$epsilon_E$*]),
    [Point / amoeba], [`(4, 0)`], [4], [$2.00 times 10^(-9)$], [$9.52 times 10^(-9)$], [$4.79 times 10^(-16)$],
    [Point / x-strip], [`(4, 0)`], [4], [$2.00 times 10^(-9)$], [$8.00 times 10^(-14)$], [$1.00 times 10^(-16)$],
    [Point / y-strip], [`(4, 0)`], [4], [$1.44 times 10^(-9)$], [$8.72 times 10^(-13)$], [$1.00 times 10^(-16)$],
    [Point / 11-strip], [`(4, 0)`], [4], [$5.55 times 10^(-10)$], [$1.44 times 10^(-9)$], [$2.10 times 10^(-16)$],
    [Line / amoeba], [`(0, 2)`], [358], [$3.48 times 10^(-10)$], [$8.45 times 10^(-7)$], [$1.32 times 10^(-10)$],
    [Line / x-strip], [`(0, 2)`], [350], [$3.48 times 10^(-10)$], [$4.22 times 10^(-7)$], [$1.32 times 10^(-10)$],
    [Line / y-strip], [`(0, 2)`], [346], [$4.41 times 10^(-10)$], [$5.80 times 10^(-7)$], [$1.28 times 10^(-10)$],
    [Line / 11-strip], [`(0, 2)`], [380], [$1.42 times 10^(-10)$], [$3.98 times 10^(-7)$], [$5.51 times 10^(-11)$],
  ),
  caption: [Comparison with the analytic GBZ radii and the fixed-energy equation for the two demonstration parameter sets. All eight cases pass the stated thresholds.],
)

Agreement means that numerical and analytic spectral membership match within `spectrum_atol` and, for nonempty results, every sampled point satisfies the closed-form GBZ constraints and fixed-energy equation within the sample tolerances. The independent membership test also checks empty results. It does not prove that every analytic branch was found in a nonempty result or that interpolation between line samples is equally accurate. An energy sweep would test additional reference energies beyond the fixed cases implemented here.


#bibliography("pyGBZ2d.bib", style: "american-physics-society")

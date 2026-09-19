#import "@preview/zhaji:0.1.0": * 

#let nt = note(
  title: "Benchmark: 2D Hatano-Nelson model",
  author: "Chenyang Wang",
  font-head: "Arial",
	font-text: "Times new roman",
	first-line-indent: 0em,
	lang: "en"
)

#show: nt.make
#show raw.where(block: true): it => block(it, fill: rgb("edf6fd"), inset: 5pt)
#show raw.where(block: false): it => highlight(it, fill: rgb("d0eafb"), top-edge: 1.2em, bottom-edge: -.5em, extent: 0.2em)

#let mathbf(x) = $bold(upright(#x))$
#let rme = $upright(e)$
#let rmi = $upright(i)$


= Benchmark: 2D Hatano-Nelson model

This example is a benchmark of the codes based on the 2D Hatano-Nelson (HN) model, where closed-form solutions are available for both amoebic GBZ and strip GBZ (see #link("https://arxiv.org/abs/2506.22743")[_arXiv_.2506.22743]). 

In this example, the following features are illustrated:
- Create a characteristic polynomial by coefficients and degrees.
- Calculate the GBZ subsets with respect to a given eigenenergy $E$, including the amoebic GBZ and the strip GBZ.
- Data readout for `LineSubset` and `PointSubset`.
- Sweep over an energy regime to get the full GBZ.

== Introduction

The Hamiltonian of the 2D Hatano-Nelson model reads,
$ H = sum_(i , j) (J_(x 1) c_(i + 1 , j)^(dagger) c_(i , j) + J_(x 2) c_(i - 1 , j)^(dagger) c_(i , j) + J_(y 1) c_(i , j + 1)^(dagger) c_(i , j) + J_(y 2) c_(i , j - 1)^(dagger) c_(i , j)) , $
where $c_(i , j)$ is the annihilation operator at the site $( i , j )$, and $J_(x 1) , J_(x 2) , J_(y 1) , J_(y 2) in CC$. In general, the coefficients can be factorized into an "Hermitian" part and a "non-Hermitian" part, i.e.,
$ J_(alpha 1) = rme^(gamma_(alpha) + rmi delta_(alpha)) J_(alpha) , quad J_(alpha 2) = rme^(- gamma_(alpha) + rmi delta_(alpha)) J_(alpha)^(*) , quad alpha = x , y $
where $gamma_(alpha) , delta_(alpha) in RR$ are the non-Hermitian factors, and $J_(alpha) in CC$ is the Hermitian factor.

In this model, the amoebic GBZ and the strip GBZs along $x$, $y$ and $[ 11 ]$ axes have #link("https://arxiv.org/abs/2506.22743")[closed forms] #cite(<wang2025generaltheorygeometrydependentnonhermitian>). The $x$-strip, $y$-strip and amoebic GBZs are identical, which read,
$ beta_(x) & = exp (gamma_(x) + rmi theta_(x)) , \
beta_(y) & = exp (gamma_(y) + rmi theta_(y)) , $
as $theta_(x) , theta_(y)$ traverses $[ 0 , 2 pi ]$. However, the $[ 11 ]$-strip GBZ is different from the above, which reads,
$ tilde(beta)_([11]) & = rme^(gamma_(x) + gamma_(y) + rmi theta_([11])) , \
tilde(beta)_(y) & = rme^(gamma_(y) + rmi theta_(y)) sqrt(lr(abs(frac(J_(x)^(*) rme^(rmi Delta_(x y) + rmi theta_([11])) + J_(y), J_(x) rme^(rmi Delta_(x y) - rmi theta_([11])) + J_(y)^(*))))) , $
where $Delta_(x y) := delta_(x) - delta_(y)$. When $sin ( Delta_(x y) ) eq.not 0$, $lr(| tilde(beta)_(y) |)$ varies with $theta_([ 11 ])$. In this case, the 2D HN model exhibits geometry-dependent bands.

The energy spectra of the strip GBZs and amoebic GBZ are illustrated in #ref(<fig:HN-model>)

#figure(
	image("Figures/HN-example-20250524.pdf"),
	caption: [ 
  Illustration of the strips and strip-GBZ spectra of the 2D HN model. *a* Illustration of the $x$-strip, $y$-strip and $[ 11 ]$-strip. *b* The spectrum of the $x$-strip and $y$-strip GBZ. The amoebic spectrum is also identical to (b). *c* The spectrum of the $[ 11 ]$-strip.   
]
)<fig:HN-model>

== Modeling tight-binding model via characteristic polynomials
`pyGBZ2d` restores a tight-binding model as a characteristic polynomial $det (E - h (beta_(1) , beta_(2)))$. Assuming that a characteristic polynomial is expanded as,
$ f (E , beta_(1) , beta_(2)) = sum_(j = 1)^(N_("terms")) a_(j) E^(l_(j)) beta_(1)^(m_(j)) beta_(2)^(n_(j)) , quad l_(j) in NN , m_(j) , n_(j) in ZZ $
then, the input format is,
```python
from pygbz2d import core
import numpy as np

coeffs = np.array([a1, a2, ..., aN])
degs = np.array([
	[l1, m1, n1],
	[l2, m2, n2],
	......,
	[lN, mN, nN],
])
poly = core.CharPoly(coeffs, degs)
```

For amoebic GBZ calculations, the sequence of the $beta$ components does not matter, since the amoebic GBZ of a given 2D model is invariant under coordinate transformations. However, strip GBZs depend on the major axis. Hence, the status of $beta_(1)$ and $beta_(2)$ are not identical for the strip-GBZ calculations. 

In `pygbz2d.sgbz` model, we assume that the input `CharPoly` is #text(c-emph)[expanded under the basis $( mathbf(a)_(1) , mathbf(a)_(2) )$, where $mathbf(a)_(1)$-axis is the major axis]. Therefore, `pygbz2d.sgbz` yields distinct results for a transformation $( beta_(1) , beta_(2) ) mapsto ( tilde(beta)_(1) , tilde(beta)_(2) )$.

For our example, we consider the basis $( mathbf(a)_(x) , mathbf(a)_(y) )$ and $( mathbf(a)_([11]) , mathbf(a)_(y) )$. For $( mathbf(a)_(x) , mathbf(a)_(y) )$, the momentum-space Hamiltonian reads,

$ h (beta_(x) , beta_(y)) = J_(x 1) beta_(x)^(- 1) + J_(x 2) beta_(x) + J_(y 1) beta_(y)^(- 1) + J_(y 2) beta_(y) , $
and the characteristic polynomial reads,
$ f (E , beta_(x) , beta_(y)) = E - J_(x 1) beta_(x)^(- 1) - J_(x 2) beta_(x) - J_(y 1) beta_(y)^(- 1) - J_(y 2) beta_(y) . $

Noted that the degrees of $beta_(x)$ and $beta_(y)$ should be permutated if you are doing $y$-strip GBZ calculation.

For $( mathbf(a)_([11]) , mathbf(a)_(y) )$, the characteristic polynomial reads,
$ f (E , beta_([11]) , beta_(y)) = E - J_(x 1) beta_([11])^(- 1) beta_(y) - J_(x 2) beta_([11]) beta_(y)^(- 1) - J_(y 1) beta_(y)^(- 1) - J_(y 2) beta_(y) . $

Detailed implementation is available in #link("./benchmark-2D-Hatano-Nelson.py")[`benchmark-2D-Hatano-Nelson.py:get_HN_charpoly`].

From the example above, you may find it disgusting to build a tight-binding model via coefficients and degrees. A recommanded approach is to use the `TightBinding` module of #link("https://github.com/Chenyang-Wang-1999/BerryPy")[#text(c-blue)[Chenyang-Wang-1999/BerryPy]], which supports easy modeling of a non-interacting tight-binding model via hopping terms. 
Besides, BerryPy.TightBinding also provides a method to calculate #text(c-emph)[the 1D GBZ] if you have #link("https://homepages.math.uic.edu/~jan/phcpy_doc_html/welcome.html")[#text(c-blue)[phcpy]] installed.
We will demonstrate the use of `BerryPy` in the application: `gain-loss-Haldane-model` and `anatomy-of-dimensional-surprise`.

== Computing GBZ subsets with specific eigenenergy $E$
Both the `sgbz` submodule and the `amoeba` submodule provide a function named `collect_GBZ_subsets` to calculate the subsets with a given eigenenergy $E$ (non-Hermitian counterpart of "equi-energy contours").
Both functions share the same prototype, which reads,
```python
def collect_GBZ_subsets(
    coeffs: np.ndarray,
    degs: np.ndarray,
    E_ref: complex,
    perc: float = None,
    **args) -> GBZResult:
    '''
    Parameters:
        coeffs: complex coefficients of the characteristic Laurent polynomial
            f(E, beta1, beta2).
        degs: (n_terms, 3) integer exponents of (E, beta1, beta2) per term.
        E_ref: reference energy to test.
        perc: progress fraction in [0, 1], printed as a percentage.
        Explicit keyword arguments only (no catch-all options dict): every
        tunable is named, defaults resolve from the home-module constants
        (see doc/constants.md), and a misspelled keyword raises TypeError
        instead of being silently ignored.

    Returns:
        GBZResult with connected subsets.
```

The class `GBZResult`, defined in `src/pygbz2d/core.py`, is a type of object that stores the results given by the solver for some eigenenergy $E_("ref")$.
The definition of `GBZResult` is,
```python
@dataclass
class GBZResult:
    """Per-E_ref GBZ computation result.

    Attributes:
        E_ref: Reference energy.
        success: Whether the computation completed without error.
        error: Error message if ``success`` is False.
        subsets: List of connected subsets (PointSubset / LineSubset).
        index: (n_0D, n_1D) counts.  ``(0, 0)`` means the energy is
               outside the GBZ.
    """

    @property
    def is_empty(self) -> bool:
        return len(self.subsets) == 0

    @property
    def is_gbz(self) -> bool:
        """True if this energy point lies on the GBZ."""
        return self.success and self.index != (0, 0)
```

Here, `GBZResult.subsets` contains all connected subsets of the contour, and `GBZResult.index` records the number of the point subsets ($n_("0D")$) and line subsets ($n_("1D")$) in `GBZResult.subsets`.

The data formats of the point and line subsets are defined as the follows,
```python
@dataclass(frozen=True)
class PointSubset:
    """0D connected subset: an isolated GBZ point.

    Represented by the (E, beta1, beta2) triplet.  Frozen (immutable)
    because it is pure data with no lazy-loaded fields.

    Attributes:
        E: Reference energy.
        beta1: Bloch factor in the principal direction.
        beta2: Representative Bloch factor in the secondary direction.
    """
    E: complex
    beta1: complex
    beta2: complex

@dataclass
class LineSubset:
    """1D connected subset: a continuous GBZ curve segment.

    Stores the curve data eagerly (no lazy fill_beta2).  One instance
    represents a single β₂ curve; n-fold degeneracies produce n LineSubsets.

    Attributes:
        E: Reference energy.
        mu1: Fixed |beta1| radius (= ln|beta1|) across the segment.
        theta1_arr: (N,) θ₁ sampling points (monotonic).
        beta2_arr: (N,) β₂ values along this single curve.
    """
    E: complex
    mu1: float
    theta1_arr: np.ndarray   # (N,)
    beta2_arr: np.ndarray    # (N,)

``` 




#bibliography("pyGBZ2d.bib", style: "american-physics-society")

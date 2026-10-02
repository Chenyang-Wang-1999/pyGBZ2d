#import "style.typ": *  // modified from @preview/zhaji:0.1.0

#let nt = note(
  title: "Geometry-dependent skin effect",
  author: "Chenyang Wang",
  font-head: "Arial",
	font-text: "Times New Roman",
	first-line-indent: 0em,
	lang: "en"
)

#show: nt.make

#let params = json("Figures/geometry-dependent-skin-effect/param.json")
#let get-param(key) = eval(params.at(key), mode: "math", scope: (
  "ii": $ii$,
  "ee": $ee$
))


= Geometry-dependent skin effect

#let fig(name, width: 100%) = image("Figures/geometry-dependent-skin-effect/" + name, width: width)
#let avec = math.bold([$a$])
#let bvec = math.bold([$b$])

This application illustrates how to identify the geometry-dependent skin effect (GDSE) via numerical calculation of the strip generalized Brillouin zone (SGBZ).
The example of interest is the gain-loss Haldane model, which is well acknowledged by its hybrid skin-topological boundary states @liGainLossInducedHybridSkinTopological2022.
It is generally recognized that the bulk bands of the gain-loss Haldane model do not exhibit non-Hermitian skin effect, as its amoebic generalized Brillouin zone (amoebic GBZ) is identical to the Brillouin zone (BZ).
However, according to its SGBZs with respect to different major axes, it is manifest that the model exhibits GDSE.

The accompany #link("./geometry-dependent-skin-effect.py")[Python script] demonstrates how to:
- Construct the characteristic polynomial via #link("https://github.com/Chenyang-Wang-1999/BerryPy")[BerryPy]
- Compute the SGBZ with a specified major axis
- Analyze the skin effect with the knowledge of the SGBZ.


== Modeling a non-Hermitian lattice with BerryPy
BerryPy is a Python package for calculating Berry phases, Berry curvature, and Chern numbers. Its `TightBinding` submodule constructs lattice Hamiltonians in a general dimension $d$ from directed hopping lists. We first introduce lattice coordinates, sublattices and coupling terms, and the momentum-space representation in $d$ dimensions. The schematics use two-dimensional examples for illustration; the final subsection specializes to the two-dimensional gain-loss Haldane model and its SGBZ calculation with pyGBZ2d.

=== Coordinate system

Consider a $d$-dimensional Bravais lattice with $d$ independent primitive vectors $avec_1,...,avec_d$. Their Cartesian components form the *columns* of the $d times d$ matrix $A=(avec_1,...,avec_d)$. A cell is labeled by an integer column vector $n=(n_1,...,n_d)^T$, and its origin is at $A n$. More generally, a real-space vector with lattice coordinates $xi$ has Cartesian coordinates $r=A xi$.

The reciprocal vectors $bvec_1,...,bvec_d$ obey $avec_i dot bvec_j=2pi delta_(i j)$ for $i,j=1,...,d$. Writing $B=(bvec_1,...,bvec_d)$ and reserving $G$ for the metric tensor of the direct lattice gives
$ G=A^T A, quad G_(i j)=avec_i dot avec_j, quad
  B=2pi A G^(-1)=2pi (A^(-1))^T, quad A^T B=2pi I_d. $
Thus $G$ contains scalar products of direct-lattice vectors, whereas the columns of $B$ are the reciprocal vectors. In BerryPy, `LatticeVec` stores $A$.

#figure(fig("direct-reciprocal-lattice.pdf", width: 80%), caption: [
  A two-dimensional illustration of primitive and reciprocal vectors for an oblique lattice. The shaded parallelograms are primitive cells. Reciprocal vectors satisfy $avec_i dot bvec_j=2pi delta_(i j)$.
]) <fig:lattice-coordinates>

A Cartesian wave vector $k=B q$ is described by fractional reciprocal coordinates $q=(q_1,...,q_d)^T$. Since $k dot avec_j=2pi q_j$, the translation phases are $theta_j=2pi q_j$ for $j=1,...,d$. The components of $q$ are defined modulo integers, and those of $theta$ modulo $2pi$. A unit $d$-dimensional hypercube in $q$, with opposite faces identified, parametrizes the BZ torus. The argument named `k` in BerryPy's `get_bulk_Hamiltonian(k)` is this fractional coordinate vector $q$, rather than the Cartesian vector $k$ used in the formulas here.

==== Coordinate-transform methods

`TightBindingModel` supplies four linear coordinate transformations. Their mathematical definitions are
$
  r=A xi, quad xi=frac(B^T,2pi)r=G^(-1) A^T r,
  quad k=B q, quad q=frac(A^T,2pi)k.
$
The corresponding interfaces are:

#table(
  columns: (auto, 1fr), inset: 6pt, stroke: .35pt + luma(75%), align: center, 
  table.header([Method], [Input and output]),
  [`lattice2cart(xi)`], [$xi arrow.r A xi$],
  [`cart2lattice(r)`], [$r arrow.r B^T r \/ 2pi$],
  [`reciprocal2cart(q)`], [$q arrow.r B q$],
  [`cart2reciprocal(k)`], [$k arrow.r A^T k \/ 2pi$],
)

A numerical input can be a single column vector, represented by a NumPy array of shape `(d,)` or `(d, 1)`, or a matrix of shape `(d, m)` whose columns are vectors. Matrix inputs are transformed column by column. A vector-valued function is transformed pointwise by composing it with the same map: for example, $k(t)=B q(t)$ and $r(t)=A xi(t)$. The current methods accept the numerical values of these functions; wrap the call to construct the transformed function, as in this straight path in fractional reciprocal coordinates:

```python
import numpy as np

q_direction = np.full(model.dim, 0.25)
q_path = lambda t: np.multiply.outer(q_direction, t)
k_path = lambda t: model.reciprocal2cart(q_path(t))
k_point = k_path(0.2)
k_columns = k_path(np.linspace(0, 1, 101))
```

Here `k_path` returns an array of shape `(d,)` for a scalar parameter or `(d, m)` for $m$ parameter samples. Passing `q_path` itself directly to `reciprocal2cart` would not perform function composition. For data stored as *rows*, transpose before and after conversion: `model.cart2lattice(points_cart.T).T`.

=== Sublattices and coupling terms

==== Basis sites and sublattice coordinates

A sublattice consists of all sites obtained by translating a chosen basis site through the Bravais lattice. A unit cell with $N$ basis sites therefore has sublattice labels $s=0,1,...,N-1$. Every translated copy of a given basis site has the same sublattice label; a site in the full lattice is specified by the pair $(n,s)$.

Let $tau_s$ be the position of basis site $s$ relative to the cell origin, expressed in the *direct-lattice basis*. Its Cartesian position is
$ r_(n,s)=A(n+tau_s). $
The integer vector $n$ is the *cell coordinate*, and the generally fractional vector $tau_s$ is the *site coordinate*; both have $d$ components. BerryPy stores these site coordinates in an $N times d$ array named `SiteCoord` (singular): row $s$ contains $tau_s^T$. Thus `SiteCoord[s]` identifies the relative position of sublattice $s$, and `SiteNum` gives $N$. The row ordering fixes the site indices used throughout the hopping lists and Hamiltonian matrices.

For Cartesian basis positions stored as rows, initialize the fractional site coordinates with

```python
model.SiteCoord = model.cart2lattice(sites_cart.T).T
```

==== Directed coupling terms

Let $t_(p s)(R)$ be the amplitude from source site $(n,s)$ to destination site $(n+R,p)$, where $R=(R_1,...,R_d)^T$ is an integer cell shift. Here $s$ and $p$ label sites, while $d$ denotes the lattice dimension. Translation invariance makes the amplitude independent of $n$. A term is *intra-cell* when $R=0$, and *inter-cell* otherwise. The physical displacement of the hopping is $A(R+tau_p-tau_s)$, but the hopping list stores only the integer shift $R$ in its cell-coordinate field.

#(nt.hint)[In the second quantization form, the hopping $t_(p s) (R)$ corresponds to the term, $ t_(p s) (R) sum_(n) c_(n+R,p)^(dagger) c_(n, s) $]

BerryPy encodes these terms as `[source, destination, amplitude]` in `InCell`, and `[source, destination, amplitude, cell_shift]` in `InterCell`, with a length-$d$ tuple for `cell_shift`. Onsite potentials have the same source and destination and belong to `InCell`. For example, the arrows in  @fig:sublattice-couplings correspond to `[0, 1, t_intra]` and `[1, 0, t_inter, (1, 0)]`. Reverse hoppings must be entered separately; they are not added automatically. Hermiticity requires $t_(s p)(-R)=overline(t_(p s)(R))$, including real onsite potentials. Violating this relation produces a non-Hermitian model.

```python
from BerryPy import TightBinding as tb

d = A.shape[1]
model = tb.TightBindingModel(
    dim=d, SiteNum=N, LatticeVec=A,
    InCell=intracell, InterCell=intercell,
    SiteCoord=site_coordinates,
)
```

#v(1em)
#figure(fig("sublattices-and-couplings.pdf", width: 90%), caption: [
  A two-dimensional, two-sublattice example with $A=I_2$. Left: cell coordinates locate cell origins, while site coordinates locate the basis sites relative to those origins; arrows distinguish intra-cell and inter-cell hopping. Right: grouping the two cells changes the site indices to $s'=0,1,2,3$, and the illustrated inter-cell hopping becomes intra-cell.
]) <fig:sublattice-couplings>
#v(1em)

==== Changing periodic units
`BerryPy.TightBindingModel` supplies re-selection of the periodic unit and basis vectors via the method:

```python
new_model = model.get_supercell(
  CoordList: list[tuple[int...]], # List of cell coordinates
  T: np.ndarray[int, 2], # Transformation matrix
)
```

Let
$ A'=A T, quad avec'_j=sum_(i=1)^d avec_i T_(i j), quad j=1,...,d, $
where $T$ is a nonsingular $d times d$ *integer* matrix. Its *columns* express the new direct-lattice vectors in the original basis. The new cell contains $L=abs(det T)$ original cells, and therefore $N L$ sites. For $L=1$, only the primitive basis changes; for $L>1$, the enlarged cell folds the bands into a smaller reciprocal cell. 

The new periodic units are specified explicitly by the cell coordinates of the included unit cells.
In `BerryPy` , `get_supercell` takes an input argument `CoordList` to record the *cell coordinates* of the merged cells.
The specified coordinates should be compatible to the transformation matrix $T$.

After setting the new cell and lattice vectors, the sublattice indices and coordinates are reassigned automatically.
Assuming the coordinates in `CoordList` are $c_l, l=0,...,L-1$, the new site index and fractional site coordinate are
$ s'=N l+s, quad tau'_(s')=T^(-1)(c_l+tau_s),
  quad l=0,...,L-1, quad s=0,...,N-1. $
In other words, BerryPy visits cells in `CoordList` order, and within each cell visits the old sites in index order. 
Row $j$ of the new `SiteCoord` contains $(tau'_j)^T$, and the full array has shape `(N * L, d)` .

In the example shown in @fig:sublattice-couplings, the supercell (b) is constructed by setting `CoordList=[(0, 0), (1, 0)]` and `T=diag([2, 1])` .
The original sites 0 and 1 in the first cell become sites 0 and 1, and those in the second cell become sites 2 and 3. 

A hopping from old sublattice $s$ in representative cell $c_l$ is reassigned by decomposing its destination as $c_l+R=c_(l')+T R'$. Its new source and destination indices are $N l+s$ and $N l'+p$, respectively. It becomes intra-cell for $R'=0$ and inter-cell otherwise, with the same amplitude. 
This relabeling retains the bonds between neighboring supercells and continues to describe an infinite periodic lattice.

=== Momentum-space Hamiltonian

The real-space eigenvalue equation associated with the directed hoppings is
$ E psi_(n,p)=sum_(s,R) t_(p s)(R) psi_(n-R,s). $
For the generalized Bloch ansatz
$ psi_(n,s)=u_s product_(j=1)^d beta_j^(n_j), $
it becomes $h(beta)u=E u$, where $beta=(beta_1,...,beta_d)$ and
$ h_(p s)(beta)=sum_R t_(p s)(R) product_(j=1)^d beta_j^(-R_j). $
A directed hopping therefore enters the *destination row and source column*, multiplied by the inverse translation factor for its cell shift.

Ordinary Bloch waves have $beta_j=ee^(ii theta_j)=ee^(2pi ii q_j)$ on the unit circle for every $j=1,...,d$. Generalized Bloch waves allow $beta_j=ee^(mu_j+ii theta_j)$ with real $mu_j$. Both `q` and `beta` below have $d$ components. The first two calls give the same sparse Hamiltonian on the BZ, while the third allows nonzero complex translation factors:

```python
H_bz = model.get_bulk_Hamiltonian(q)
H_bz_equivalent = model.get_bulk_Hamiltonian_complex(
    np.exp(2j * np.pi * np.asarray(q)))
H_nonbloch = model.get_bulk_Hamiltonian_complex(beta)
```

==== Cell and site phase conventions

`get_bulk_Hamiltonian_complex` uses only the integer cell coordinates in its Bloch factors. The entries of `SiteCoord` do not enter this Hamiltonian: the phase of a hopping with cell shift $R$ is $ee^(-ii k dot (A R))$, rather than $ee^(-ii k dot (A(R+tau_p-tau_s)))$. These two phase conventions are related by a change of basis.

To see this, write the same physical wavefunction in the two conventions:
$ psi_(n,s)=ee^(ii k dot (A n)) u_s^"cell"
  =ee^(ii k dot (A(n+tau_s))) u_s^"site". $
Define the diagonal matrix
$ D_(s s)(k)=ee^(ii k dot (A tau_s))=ee^(2pi ii q^T tau_s). $
Then
$ u^"cell"=D(k)u^"site", quad
  h_"site"(k)=D(k)^(-1) h_"cell"(k)D(k), $
and the matrix elements in the site convention are
$ (h_"site")_(p s)(k)=sum_R t_(p s)(R)
  ee^(-ii k dot (A(R+tau_p-tau_s))). $
Thus each site amplitude undergoes a phase rotation determined by its own intracell position. Because `SiteCoord[s]` is fractional, the phase is $2pi q^T tau_s$, or equivalently $k dot (A tau_s)$ with Cartesian $k$. At real momentum $D$ is unitary, and for every fixed momentum the two Hamiltonians have the same eigenvalues and characteristic polynomial by similarity.

The equivalence also extends to non-Bloch factors on a chosen local branch of $log beta_j$:
$ D_(s s)(beta)=exp(sum_(j=1)^d (tau_s)_j log beta_j), quad
  h_"site"(beta)=D(beta)^(-1)h_"cell"(beta)D(beta). $
Here $D$ is generally nonunitary but remains invertible for nonzero $beta_j$. The cell convention used by BerryPy keeps the Hamiltonian Laurent-polynomial in the translation factors.

==== Opening selected boundary directions

In `get_bulk_Hamiltonian_complex`, any subset of the $d$ directions can be opened. With mathematical axes numbered $j=1,...,d$, setting the Python entry `beta[j - 1]` to `None` removes every intercell hopping with $R_j != 0$. Hoppings within the current cell are retained.

Let $cal(O)$ contain the open-axis indices and $cal(P)$ the remaining indices. Denote by $cal(S)$ the integer shifts with $R_j=0$ for every $j in cal(O)$. The resulting Hamiltonian is
$ h_(p s)^"open"(beta)=sum_(R in cal(S)) t_(p s)(R)
  product_(j in cal(P)) beta_j^(-R_j). $
When all $d$ directions are open, this reduces to $h_(p s)^"open"=t_(p s)(0)$. In the two-dimensional illustration,
$ h_(p s)(beta_1,"None")=sum_(R:R_2=0)t_(p s)(R)beta_1^(-R_1),
  quad h_(p s)("None","None")=t_(p s)(0). $
Here $t_(p s)(0)$ denotes all intracell contributions in the *current* model, including the bonds between old cells that have been grouped into a supercell. `None` removes bonds; it is not a numerical value of $beta$ or the limit $beta arrow.r 0$.

#figure(fig("boundary-conditions.pdf"), caption: [
  A two-dimensional illustration using the same $4 times 3$ supercell with both directions periodic, with direction 2 open, and with both directions open. Blue bonds connect neighboring supercells; red crosses mark the bonds discarded by `None`. Internal bonds and site numbers are unchanged. Unit-modulus values of the remaining $beta_j$ give ordinary Bloch boundary conditions in the retained directions.
]) <fig:boundary-conditions>

To obtain a sample with chosen widths along the open directions, first expand the cell to contain those sites, then apply `None` to the expanded model. For a two-dimensional, one-site square-lattice model, the example in @fig:boundary-conditions is

```python
cells = [(i, j) for i in range(4) for j in range(3)]
supercell = model.get_supercell(cells, np.diag([4, 3]))
beta1, beta2 = np.exp(2j * np.pi * np.array([0.2, 0.3]))
H_periodic = supercell.get_bulk_Hamiltonian_complex((beta1, beta2))
H_strip = supercell.get_bulk_Hamiltonian_complex((beta1, None))
H_open = supercell.get_bulk_Hamiltonian_complex([None] * supercell.dim)
```

The remaining factors refer to translations of the *supercell*. In any dimension, a length-$d$ list of `None` opens all directions. If boundary conditions are specified using fractional momenta, a `None` component has the same intended meaning. The current `get_bulk_Hamiltonian(k)` wrapper assumes numerical components, so preserve `None` explicitly when converting to the complex-factor interface. The conversion below works for a list `q` of any length $d$; the displayed values continue the two-dimensional example:

```python
q = (0.2, None)  # fractional momentum in direction 1; open direction 2
beta = [None if qj is None else np.exp(2j * np.pi * qj) for qj in q]
H_strip = supercell.get_bulk_Hamiltonian_complex(beta)
```

For a fully periodic finite-range model, the characteristic equation
$ f(E,beta_1,...,beta_d)=det(E I_N-h(beta))=0 $
is a Laurent-polynomial equation. `model.get_characteristic_polynomial_data()` returns its coefficients and exponent rows in the order $(E,beta_1,...,beta_d)$, with $d+1$ exponents per row.

Under $A'=A T$, a generalized Bloch wave has translation factors
$ beta'_j=product_(i=1)^d beta_i^(T_(i j)), quad j=1,...,d, quad
  mu'=T^T mu, quad theta'=T^T theta quad (mod 2pi). $
The unit-modulus condition defining the BZ is preserved by a nonsingular $T$ in any dimension. For the two-dimensional SGBZ calculation in this application, construct a cell whose first vector $avec'_1$ lies along the chosen boundary and extract the characteristic polynomial *with the intercell hoppings retained*. The ordering of the two transformed axes sets the major and minor directions for pyGBZ2d. The solver then determines the generalized Bloch factors; the explicit `None` construction serves finite-matrix calculations.

#(nt.hint)[The real-space Hamiltonian in the parallelogram geometry can be easily constructed by combining `get_supercell` and `get_bulk_Hamiltonian((None, None))` .
However, this method is inefficient
]

=== Modeling the gain-loss Haldane model

We now specialize to $d=2$ and apply these conventions to a honeycomb lattice with two orbitals, labeled $A$ and $B$ (indices 0 and 1). The model has real nearest-neighbor hopping $t_1$, complex next-nearest-neighbor hopping with magnitude $t_2$ and opposite circulating phases on the two sublattices, and onsite potentials $M$ and $-M$. Taking $M$ to be purely imaginary introduces balanced gain and loss. The parameter choice used below has Hermitian hopping terms; its non-Hermiticity comes from the onsite potentials.

#figure(fig("geometries.pdf"), caption: [
  Two open geometries built from the same hopping model. Red and blue sites carry gain and loss, respectively. Arrows denote the lattice directions used by the script; their positive senses need not coincide with positive Cartesian axes. The rhombus has sides parallel to $avec_1$ and $avec_2$; the rectangle has sides parallel to $x$ and $y$. The drawings show small lattices; the OBC data below use 12,800 sites in each geometry.
]) <fig:geometry>


Use a two-site primitive cell with lattice vectors and Cartesian site positions
$
  A = (avec_1, avec_2) = mat(-1/2, -1/2; -sqrt(3)/2, sqrt(3)/2), quad
  r_A = (0, 1/(2 sqrt(3))), quad r_B = (0, -1/(2 sqrt(3))).
$
The parameters follow the supplied Haldane scripts:
$ t_1=#get-param("t1"), quad t_2=#get-param("t2"), quad phi=#get-param("phi"), quad M=#get-param("M"), quad gamma=#get-param("gamma") $
Here the gain-loss strength is $"Im" M=0.5$. The additional parameter $gamma$ in the script is a common phase multiplying both directions of the next-nearest-neighbor hopping; it is set to zero in this case.

The following is the $gamma=0$ construction used by `build_model()`. Each nearest-neighbor bond and its reverse are entered explicitly. For the next-nearest-neighbor bonds, `sign` distinguishes the two sublattices and the reversed cell shift carries the conjugate hopping phase.

```python
t1, t2, phi, M = 1.0, 0.5, np.pi / 3, 0.5j
A = np.array([[-0.5, -0.5],
              [-np.sqrt(3) / 2, np.sqrt(3) / 2]])
intracell = [[0, 0, M], [1, 1, -M],
             [1, 0, t1], [0, 1, t1]]
intercell = [[1, 0, t1, (0, -1)], [1, 0, t1, (1, 0)],
             [0, 1, t1, (0, 1)], [0, 1, t1, (-1, 0)]]
for s, sign in ((0, 1), (1, -1)):
    for R in ((-1, 0), (0, -1), (1, 1)):
        intercell.append([s, s, t2 * np.exp(1j * sign * phi), R])
        intercell.append([
            s, s, t2 * np.exp(-1j * sign * phi),
            tuple(-r for r in R),
        ])
model = tb.TightBindingModel(2, 2, A, intracell, intercell)
sites_cart = np.array([[0, 1 / (2 * np.sqrt(3))],
                       [0, -1 / (2 * np.sqrt(3))]])
model.SiteCoord = model.cart2lattice(sites_cart.T).T
```

Using the hopping-to-matrix rule above, define
$ Q=beta_1+beta_2+(beta_1 beta_2)^(-1), quad tilde(Q)=beta_1^(-1)+beta_2^(-1)+beta_1 beta_2. $
The resulting two-band Hamiltonian is
$ h(beta_1,beta_2)=mat(
  M+t_2(ee^(ii phi)Q+ee^(-ii phi)tilde(Q)), t_1(1+beta_2+beta_1^(-1));
  t_1(1+beta_2^(-1)+beta_1), -M+t_2(ee^(-ii phi)Q+ee^(ii phi)tilde(Q))
). $
On the BZ, $tilde(Q)=overline(Q)$. Away from the unit torus, $tilde(Q)$ remains the explicitly defined Laurent expression above; it is not obtained by complex-conjugating $Q$. This distinction preserves the analytic continuation needed for the SGBZ calculation. The script samples the BZ using `get_bulk_Hamiltonian_complex(np.exp(1j * theta))`, where `theta` contains phases in radians.

We next choose the periodic cell for each strip orientation. The primitive directions $avec_1$ and $avec_2$ and the horizontal direction $x$ are parallel to zigzag boundaries; the vertical direction $y$ is parallel to armchair boundaries. The first vector in each row below is the strip's major axis. Here $x$ and $y$ name the Cartesian axis directions; the positive senses of the chosen lattice vectors follow the signs in the table.

#table(
  columns: (auto, 1fr, 1fr, auto), inset: 6pt, stroke: .35pt + luma(75%),
  table.header([Name], [Major axis $avec'_1$], [Minor axis $avec'_2$], [Sites/cell]),
  [`a1`], [$avec_1$], [$avec_2$], [2],
  [`a2`], [$avec_2$], [$avec_1$], [2],
  [`x`], [$avec_1+avec_2$], [$avec_1-avec_2$], [4],
  [`y`], [$avec_1-avec_2$], [$avec_1+avec_2$], [4],
)

```python
cells = [(0, 0), (1, 0)]
T_x = np.array([[1, 1], [1, -1]])
T_y = T_x[:, [1, 0]]
model_x = model.get_supercell(cells, T_x)
model_y = model.get_supercell(cells, T_y)
coeffs, degs = polynomial_data(model_x)
```

Both Cartesian cells contain two primitive cells because $abs(det T_x)=abs(det T_y)=2$. For `a1`, `direction_model()` uses the original two-site model. For `a2`, it uses the single representative `(0, 0)` and the matrix `[[0, 1], [1, 0]]` to exchange the primitive axes. The four-site `x` and `y` models exchange their major and minor axes in the same way. Their phase plots consequently use different coordinates and folded bands. Extracting the characteristic polynomial after each transformation gives the SGBZ solver the hopping data for the intended boundary orientation.



== Scanning and displaying the SGBZ

For each direction, the script obtains the Laurent coefficients of
$ f(E,beta_1,beta_2)=det(E I-h(beta_1,beta_2)) $
from BerryPy. Each row of `degs` contains the exponents of $(E,beta_1,beta_2)$ in that order. It then scans a rectangular mesh in the complex energy plane:

```python
from pygbz2d.sgbz import collect_GBZ_subsets

coeffs, degs = polynomial_data(direction_model("y"))
result = collect_GBZ_subsets(coeffs, degs, 1.5 + 0j)
if not result.success:
    print(result.error)             # unresolved, not outside the spectrum
elif result.is_gbz:
    print(result.index, result.subsets)
else:
    print("No SGBZ subset at this energy")
```

The solver locates a zero of the average major-axis winding and imposes the transverse equal-modulus root condition. A vanishing winding plateau is excluded by the solver's plateau check. This example uses the package defaults; it does not impose $mu_1=0$, set transverse radii to one, or remove discrepant samples.

#figure(fig("sgbz-comparison.svg"), caption: [
  Directional comparison using the existing $201 times 201$ energy scans, with $"Re" E in [-3.1,4.6]$ and $"Im" E in [-0.51,0.51]$. Top: sampled SGBZ spectra over the BZ spectrum in gray. Bottom: phases $(theta_1,theta_2)$ of the computed subsets, colored by $max_j abs(ln abs(beta_j))$. Every stored line-subset row is included. White regions indicate missing samples in this finite energy scan, not established holes in the BZ or SGBZ. Each panel uses its own lattice basis.
]) <fig:sgbz>

`PointSubset` supplies one pair $(beta_1,beta_2)$. For a `LineSubset`, the script uses *all* entries of `theta1_arr` and `beta2_arr`, reconstructing $beta_1=exp(mu_1+ii theta_1)$. This matters for the real-energy continuum subsets in the $y$ scan. Failed solves and unfinished checkpoints remain distinguishable from successfully classified spectral exterior points.

#pagebreak()
== Reading the skin-effect diagnostic

Write $beta_j=exp(mu_j+ii theta_j)$. A nonzero $mu_j=ln abs(beta_j)$ describes exponential growth or decay per cell in direction $avec'_j$. The BZ is $mu_1=mu_2=0$. Checking only the major-axis radius would miss the effect here: the $y$-SGBZ has $mu_1$ close to zero while its transverse $mu_2$ is generically nonzero. For this strip, $avec'_2$ points along $x$, so the localization is transverse to the armchair boundary.

#figure(fig("sgbz-radii.svg"), caption: [
  Logarithmic radii of every stored subset sample, plotted against its real energy. The small scale of the major-axis panel reflects numerical deviations near zero. The broad nonzero transverse radii in the armchair case are the main diagnostic; the sample distribution in these plots is not a density of states.
]) <fig:radii>

The following statistics are generated by `sweep_report` and stored in `sgbz-summary.json`. Each scan contains 40,401 energies, with no entries marked failed or unfinished.

#table(
  columns: (auto, auto, 1fr, 1fr, 1fr), inset: 5pt, stroke: .35pt + luma(75%),
  table.header([Direction], [In spectrum], [Max. $abs(mu_1)$], [Max. $abs(mu_2)$], [99th % $abs(mu_2)$]),
  [$avec_1$], [9,991], [$2.2 times 10^(-16)$], [$9.15 times 10^(-10)$], [$2.62 times 10^(-11)$],
  [$avec_2$], [9,991], [$2.2 times 10^(-16)$], [$9.15 times 10^(-10)$], [$2.62 times 10^(-11)$],
  [$x$], [9,991], [$1.46 times 10^(-5)$], [$1.15 times 10^(-4)$], [$4.58 times 10^(-11)$],
  [$y$], [6,656], [$1.00 times 10^(-4)$], [$0.590$], [$0.505$],
)

These samples support a BZ description of the bulk for the three zigzag directions and a non-Bloch description for the armchair direction. The latter is not merely an energy-spectrum comparison: the transverse Bloch factors themselves have nonunit modulus. A finite scan does not establish an exact equality of continuous manifolds, nor does it establish the behavior of every possible polygon.

=== Numerical deviations retained in the figures

At $E=1.674+0.2397 ii$, the legacy $x$ scan has six point subsets with $mu_1 approx 1.46 times 10^(-5)$ and maximum $abs(mu_2) approx 1.15 times 10^(-4)$. A fresh call to the current solver reproduces this deviation. A fixed-$mu_1=0$ winding evaluation is also slightly nonzero (about $-1.55 times 10^(-5)$), so this is not simply a plotting or pickle-conversion problem. Its full numerical cause remains unresolved. The figure and maximum statistic retain it; the 99th percentile shows how localized it is in this scan. The small $y$ major-axis deviations are likewise retained.

#pagebreak()
== Constructing finite open geometries

The rhombus repeats the two-site primitive cell over $0 <= n_1<N_1$, $0 <= n_2<N_2$. Its sides are parallel to $avec_1$ and $avec_2$. The rectangle repeats the four-site `x` cell over the same index ranges. Its vertical sides are armchair boundaries, and its horizontal sides are zigzag boundaries. The choice is guided by the directional SGBZ comparison, then checked with finite OBC eigenstates.

`finite_hamiltonian` reads BerryPy's intracell and intercell hopping lists. It retains a bond only when both endpoint cells belong to the finite sample. The right eigenvectors are obtained using the general complex eigensolver `scipy.linalg.eig`; the matrix is non-Hermitian, so `eigh` would be inappropriate.

#figure(fig("obc-comparison.svg"), caption: [
  Existing OBC data for an $80 times 80$ primitive-cell rhombus and an $80 times 40$ four-site-cell rectangle, both containing 12,800 sites. Top: finite OBC spectra and the BZ reference. Gap states are included in these spectra. Bottom: mean normalized right-state intensity in the shaded window $-2.5 <= "Re" E <= -1.5$, containing 5,178 and 5,231 states, respectively. Both maps share a logarithmic color scale; a uniform density would have value one. The rectangle exhibits pronounced accumulation near its armchair sides, whereas the rhombus retains substantial bulk weight. Finite-size modulations and corner enhancement remain visible in the rhombus.
]) <fig:obc>

The finite spectra need not coincide point by point with the infinite-strip spectra: the order of limits differs, and finite OBC spectra also contain boundary states. Their role here is to check the geometric construction and the spatial localization predicted by the directional comparison.

#pagebreak()
== State density and data provenance

Every right eigenvector is normalized separately before averaging. For the selected set $cal(W)$, the plotted quantity is
$ rho(r)=frac(N, abs(cal(W))) sum_(n in cal(W)) frac(abs(psi_(R,n)(r))^2, sum_(r') abs(psi_(R,n)(r'))^2), quad
  frac(1,N) sum_r rho(r)=1. $
The energy window lies inside the lower bulk band and avoids the central topological gap. It is a reproducible energy selection, not a rigorous topological classification of each finite-system state. Averaging many such states reduces the risk of mistaking a single topological corner mode for extensive bulk skin localization. The compact archives also retain the all-state mean, inverse participation ratios, and one deterministically selected state (nearest to $-2+0.2 ii$ within the same window).

=== Sources of the figures

- *SGBZ:* the four supplied files `data/Haldane-gain-loss-{a1,a2,x,y}-SGBZ.pkl`. The script checks the energy-grid ordering and model parameters when reading them. The two-site scans have 26 polynomial terms; the four-site scans have 64.
- *OBC rhombus:* `paper-Haldane-gain-loss-OBC-80-80.pkl` from `myoffice:~/654/research-data/2D_skin_effect/phcpy-free-algorithm/data`.
- *OBC rectangle:* `paper-Haldane-gain-loss-OBC-square-80-40.pkl` from the same directory. The legacy name “square” refers to the Cartesian cell construction; the actual finite sample need not be a geometric square.
- *Compact extracts:* `obc-rhombus.npz` and `obc-rectangle.npz` beside the figures. Each retains all eigenvalues and site coordinates together with normalized density summaries. The original eigenvector arrays are not duplicated. Source paths, file sizes, modification times, and extraction times are stored in each archive's `metadata` field.

The legacy OBC tuples contain no model-parameter metadata. Their parameter assignment follows the companion Haldane script supplied with those calculations. The original caches are preserved; the figures represent reused research data, not fresh large-scale calculations in this example.

=== Validation performed for this example

The `check` command compares the BerryPy characteristic polynomial with $det(E I-h)$ at eight deterministic complex arguments per direction. Maximum relative errors are below $3.3 times 10^(-15)$. Independently, direct finite-lattice assembly is compared with BerryPy's two-stage supercell construction followed by `get_bulk_Hamiltonian_complex((None, None))`: all four matrices agree exactly, with coordinate differences below $4.5 times 10^(-16)$.

The small demonstration samples six energies per direction, including spectral interior, exterior, and a real-energy continuum case, and diagonalizes two 128-site OBC samples. It tests the complete pipeline; those small lattices are not used to claim thermodynamic convergence. The detailed validation record, including unexpected results, accompanies this case in `log/2026-09-25-geometry-dependent-skin-effect.md`.

#pagebreak()
== Running and extending the example

Run the commands below from the repository root. Multiline examples use POSIX shell continuation (`\`); in PowerShell, enter the command on one line or use its backtick continuation. The script resolves default directories relative to its own location. Install this checkout with `python -m pip install -e .`, and use an environment containing BerryPy, SymPy, NumPy, SciPy, and Matplotlib. These application dependencies are separate from the lightweight `pygbz2d` package dependencies.

=== Small end-to-end calculation

```sh
python application/geometry-dependent-skin-effect.py demo --workers 2
python application/geometry-dependent-skin-effect.py check
```

The demo writes to `application/data/geometry-dependent-skin-effect/demo`. Its grid is $"Re" E in {-2,1.5,5}$ and $"Im" E in {0,0.3}$. Intermediate results are saved after each energy. A repeat invocation reuses completed SGBZ results; OBC diagonalization and figures are regenerated.

=== Larger SGBZ scans on a server

```sh
python application/geometry-dependent-skin-effect.py sweep \
  --directions a1 a2 x y --n-real 201 --n-imag 201 --workers 12
```

The default sweep directory is `application/data/geometry-dependent-skin-effect/sweep`. The default energy bounds match the old scans. Reuse exactly the same command to resume missing entries; add `--retry-failed` to recompute failed entries. A checkpoint with incompatible direction, grid, polynomial, or parameters is rejected. Choose another `--output-dir` when changing the run. The default `sweep` grid, if sizes are omitted, is only $7 times 5$.

=== Finite OBC calculations

```sh
python application/geometry-dependent-skin-effect.py obc \
  --geometry rhombus --nx 80 --ny 80 \
  --output application/data/geometry-dependent-skin-effect/rhombus.npz
python application/geometry-dependent-skin-effect.py obc \
  --geometry rectangle --nx 80 --ny 40 \
  --output application/data/geometry-dependent-skin-effect/rectangle.npz
```

Add `--save-vectors` to retain the full eigensystem in a neighboring `.pkl` file. Without it, only the compact summary is saved. Dense diagonalization scales cubically in the number of sites, and one $12,800 times 12,800$ complex matrix alone occupies about 2.44 GiB; the eigensolver needs additional workspace. Start with `--nx 8 --ny 8` when checking a new environment.

#pagebreak()
== Reusing results and rebuilding this note

=== Extracting old OBC data where it is stored

The `summarize-obc` command reads the tuple `(eigenvalues, right_eigenvectors, coordinates)` produced by the old script. It processes eigenvector columns in batches to avoid allocating another full square probability matrix. This extraction path needs NumPy but does not import BerryPy or `pygbz2d`.

```sh
python application/geometry-dependent-skin-effect.py summarize-obc \
  /path/to/paper-Haldane-gain-loss-OBC-80-80.pkl \
  --output obc-rhombus.npz
```

Run extraction on the data server and copy the small `.npz` result for plotting. For streaming over SSH, `--output -` writes a binary NPZ archive to standard output. Use a binary-safe client when saving that stream. Pickle inputs must come from trusted calculations.

=== Recreating the figures from the supplied caches

```sh
python application/geometry-dependent-skin-effect.py plot --obc \
  application/Figures/geometry-dependent-skin-effect/obc-rhombus.npz \
  application/Figures/geometry-dependent-skin-effect/obc-rectangle.npz
```

This reads the SGBZ files from the repository's `data` directory and writes SVG, PNG, and a JSON numerical summary to `application/Figures/geometry-dependent-skin-effect`. Add `--sgbz-dir application/data/geometry-dependent-skin-effect/sweep` to visualize a newly computed scan. `--directions y` restricts the plot to one direction; `--nk` controls the BZ reference sampling. Plotting performs no SGBZ root solving or OBC diagonalization.

```sh
typst compile application/geometry-dependent-skin-effect.typ
```

The compiled note uses the checked-in figures. For a convergence study, increase the energy-grid resolution and the OBC sizes separately, retain the same density window across sizes, and inspect both logarithmic radii and real-space distributions. Small isolated radius deviations should remain visible in the numerical record until their cause is understood.

== References
#bibliography("pyGBZ2d.bib", title: none, style: "american-physics-society")

Online versions: #link("https://doi.org/10.1103/PhysRevLett.128.223903")[Phys. Rev. Lett. 128, 223903 (2022)]; #link("https://arxiv.org/abs/2506.22743v3")[arXiv:2506.22743v3].



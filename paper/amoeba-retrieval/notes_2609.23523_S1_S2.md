# Extraction: arXiv:2609.23523v1 — "Universal Generalized Brillouin Zone Theory I: Review of the Spectral Approach"

**Assigned reading region = Abstract + Section I (incl. I.1, I.2) + Section II.**

- Authors: Zeqi Xu, Jiangping Hu, Zhesen Yang
- Affiliations: (1) Dept. of Physics, Xiamen University, Xiamen 361005, Fujian, China; (2) Beijing National Laboratory for Condensed Matter Physics and Institute of Physics, Chinese Academy of Sciences, Beijing 100190, China; (3) New Cornerstone Science Laboratory, Beijing 100190, China; (4) Asia Pacific Center for Theoretical Physics, Pohang 37673, Korea
- Corresponding authors / e-mails seen in the HTML author block: `jphu@iphy.ac.cn`, `yangzs@xmu.edu.cn`
- arXiv metadata (from https://arxiv.org/abs/2609.23523v1): `arXiv:2609.23523v1 [cond-mat.mes-hall]`, submitted **Sun, 20 Sep 2026 10:18:51 UTC**, comments field: **"16 pages, 10 figures"**, DOI `10.48550/arXiv.2609.23523`, license "arXiv.org perpetual, non-exclusive license", source size 14,721 KB, HTML is "HTML (experimental)".
- Paper-internal date line: "September 20, 2026".
- HTML source used: https://arxiv.org/html/2609.23523v1 (fetched twice, plus `#S2` fragment; identical content and identical truncation point each time).

---

## 0. FETCH / TRUNCATION REPORT (honest statement of what I could and could not read)

- `web_fetch` of `https://arxiv.org/html/2609.23523v1` **always truncates at the same point**: inside the *unnumbered* real-space OBC eigenvalue equation of Section II (the large block-matrix equation, in the middle of the RHS column vector's LaTeX, right after the fragment `... \hbox to0pt{\lxSVG@begingroup@{\_sc`). The tool message is explicit: "*(Content truncated. Fetch a more specific URL or section for the full text.)*"
- Re-fetching with the `#S2` fragment returned **byte-identical** content and the **same** truncation point (the fetch is not fragment-aware).
- `https://www.arxiv.org/html/2609.23523v1` and `https://ar5iv.labs.arxiv.org/html/2609.23523` both failed: "cross-origin redirect to https://arxiv.org is not followed automatically".
- `https://arxiv.org/pdf/2609.23523v1` failed: "unsupported content type application/pdf" (no PDF text extraction available to me).
- Therefore, **exactly readable**: Abstract, Section I (all of it, including I.1 and I.2), and Section II **only up to and including the introduction of the OBC real-space eigenvalue equation** (i.e. Eqs. (10), (11), (12) and the sentence initiating the block matrix equation).
- **NOT readable / not present in fetched region**: the remainder of Section II (everything after the block-matrix equation), and Sections III–VIII and the References list (only the hyperlinked table of contents for them is visible).
- Sections III–VIII titles ARE visible via the HTML navigation/TOC (see §7 below) — titles only, no body text.

MathJax markup in the fetch is doubled/escaped (e.g. `|\\beta\_{p}|=\\beta\_{p+1}`); all equations below are reproduced **cleaned up** (un-doubled), which is a pure notational cleanup and adds no content. All *prose* quoted verbatim is given exactly as it appeared, with the same doubled-macro ugliness removed only inside math.

---

## 1. Abstract (verbatim, math cleaned)

> "This series of papers aims to establish a universal generalized Brillouin zone (GBZ) theory for higher-dimensional non-Hermitian systems. As the starting point of this series, we emphasize a fundamental question: while the conventional one-dimensional (1D) GBZ condition, |\beta_p|=|\beta_{p+1}|, is well known to fail in two dimensions, how does this breakdown actually occur as a system gradually crosses over from 1D to 2D? Investigating this question reveals that the existing 1D GBZ theory itself remains incomplete. In this first paper, we therefore systematically review the 1D spectral approach and clarify where the underlying difficulties lie. Our work identifies the key challenges that motivate the wavefunction approach developed in Paper II."

Note: the abstract writes the conventional 1D GBZ condition **without the prefactor indices** (`|\beta_p|=|\beta_{p+1}|`), whereas in the body of Sec. I it is written with band/enlarged-count indices as `|\beta_{np}(E)|=|\beta_{np+1}(E)|` (Eq. 4). Both forms appear verbatim; I record the discrepancy rather than harmonize it.

---

## 2. SECTION I — Introduction

### 2.1 Untitled preamble of Sec. I (computational framing)

Framing of the computational object: the OBC spectrum/wavefunctions come from diagonalizing H_OBC.

> "The energy spectrum and wavefunctions provide the starting point for understanding the physical properties of a lattice model. For a finite system with open boundary conditions (OBCs), they can be obtained by diagonalizing its OBC Hamiltonian H_{\mathrm{OBC}}. With L unit cells and n internal degrees of freedom per cell, this is an nL\times nL matrix. Therefore, direct diagonalization becomes impractical for macroscopically large systems, making analytical methods essential."

> "For periodic Hermitian lattices, Bloch's theorem makes this problem tractable. Under periodic boundary conditions (PBCs), one need only diagonalize the n\times n Bloch Hamiltonian H(\bm{k}) at each crystal momentum \bm{k}."

> "To apply this description to an open system, one assumes that changing the boundary conditions from PBC to OBC leaves the bulk spectrum unchanged in the thermodynamic limit [1]. Although this is generally true for Hermitian systems, it can fail in non-Hermitian systems. A prominent example is the non-Hermitian skin effect (NHSE), in which a macroscopic number of eigenstates become exponentially localized near the boundaries."

> "In such systems, the PBC and OBC spectra can differ substantially in the complex energy plane."

> "Therefore, the NHSE calls for an extension of the Bloch description. To accommodate exponential localization, the Brillouin zone (BZ) is extended from real to complex momenta [37]. In one dimension (1D), we write the corresponding translation factor as"

**Eq. (1)** (numbered):
```
\beta \equiv e^{ik} e^{\mu},
```

> "where k and \mu are real. The phase k describes spatial oscillations, while \mu=\ln|\beta| sets the rate of exponential growth or decay."

> "In the conventional BZ under PBCs, we have \mu=0, so that \beta traces the unit circle. However, under OBCs, \mu generally depends on k and need not vanish. This deformation naturally leads to the notion of the generalized Brillouin zone (GBZ), defined as the set of complex translation factors that reproduces the bulk OBC spectrum in the thermodynamic limit [37, 41, 35]. If the GBZ radius is a single-valued function of k, the GBZ of the m-th band can be written as"

**Eq. (2)** (numbered):
```
\beta_{m,\rm GBZ}(k) \equiv e^{ik} e^{\mu_{m,\rm GBZ}(k)},
```

> "where m=1,\ldots,n labels the bands. The bulk OBC spectrum is then obtained by evaluating the corresponding dispersion relation along each GBZ:"

**Eq. (3)** (numbered):
```
\sigma_{\mathrm{OBC,bulk}}=\bigcup_{m}\left\{E_{m}\big(\beta_{m,\rm GBZ}(k)\big)\,\Big|\,k\in[0,2\pi)\right\}.
```

> "Here E_{m}(\beta) denotes the dispersion branch associated with the m-th GBZ. **The central question is therefore: what determines \mu_{m,\rm GBZ}(k)?**"

> "The standard approach starts from the characteristic equation (ChE), \det[E\mathbbm{I}_{n}-H(\beta)]=0. We order its roots by increasing modulus:"

**Unnumbered display equation** (appears in a table row with an empty equation-number cell, i.e. no number visible in the fetched HTML):
```
|\beta_{1}(E)| \leq |\beta_{2}(E)| \leq \cdots \leq |\beta_{n_{s}}(E)|,
```

> "where n_{s} is the number of roots, counted with multiplicity. For an n-band model with full-rank longest-range hopping matrices, n_{s}=n(p+q), where p and q are the maximum hopping ranges to the right and left, respectively. For generic models satisfying the full-rank condition discussed in Sec. II, the bulk OBC spectrum in the thermodynamic limit is characterized by the conventional GBZ condition [37, 41, 35]:"

**Eq. (4)** (numbered) — the central object of the whole series:
```
|\beta_{np}(E)|=|\beta_{np+1}(E)|.
```

> "This condition describes the thermodynamic bulk spectrum and does not by itself determine finite-size corrections or the quantization of topological edge states. Even in this limit, however, exceptions are known [13, 38]. **Understanding when and why the conventional condition fails is a central goal of this series.**"

Additional quotes worth preserving:

> "The central question is therefore: what determines \mu_{m,\rm GBZ}(k)?"  (already above)

> "The 1D GBZ theory is well established for conventional cases [37, 41, 35, 38, 45, 24]. However, extending it to two dimensions (2D) and beyond remains a central problem, with several approaches addressing different aspects [14, 46, 33, 11, 42, 43, 12, 7, 30, 10, 29, 9, 44, 31, 32, 25, 28]. This series aims to develop a unified GBZ theory for higher-dimensional systems and clarify its relation to other approaches. A detailed review of existing 2D theories will be given in subsequent papers." (opening of §I.1)

### 2.2 Section I.1 — Motivating example (THE core numerical setup, quoted in full detail)

> "Although there is no well-accepted 2D GBZ theory, it is clear that the 1D GBZ condition in Eq. (4) cannot be generalized to higher dimensions. **Therefore, our motivating question is: how does the 1D GBZ condition in Eq. (4) break down as the system crosses over from 1D to 2D, and how can this breakdown be characterized?**"

> "As a concrete example, we consider a 2D single-band tight-binding model with the non-Bloch Hamiltonian"

**Eq. (5)** (numbered) — the model Hamiltonian:
```
H(\beta_{x},\beta_{y}) = t_{x}\beta_{x} + t_{-x}/\beta_{x} + t_{y}\beta_{x}\beta_{y} + t_{-y}/(\beta_{x}\beta_{y}).
```

> "Here \beta_{x} and \beta_{y} are the translation factors in the two lattice directions, and t_{\pm x} and t_{\pm y} are the hopping amplitudes. We impose OBCs in both directions on a rectangular lattice with L_{x} sites along x and n_{y} sites along y. We group the n_{y} sites at each x coordinate into one supercell, reducing the problem to an effective 1D model. Increasing n_{y} can be regarded as a crossover from 1D to 2D."

> "Therefore, replacing translation along x by the complex factor \beta_{x}, the transverse chain has the hopping amplitudes"

**Eq. (6)** (numbered; one equation, three lines — the number appears only on the first line in the HTML table):
```
T_{y,0}(\beta_{x})  &= t_{x}\beta_{x} + t_{-x}/\beta_{x}, \\
T_{y,+1}(\beta_{x}) &= t_{y}\beta_{x}, \\
T_{y,-1}(\beta_{x}) &= t_{-y}/\beta_{x},
```

> "where T_{y,0} is the onsite energy and T_{y,\pm1} multiply (\beta_{y})^{\pm1}. The resulting n_{y}\times n_{y} Hamiltonian reads"

**Eq. (7)** (numbered) — the supercell (transverse) Hamiltonian:
```
H_{n_{y}}(\beta_{x}) =
\begin{pmatrix}
T_{y,0} & T_{y,+1} & & \\
T_{y,-1} & T_{y,0} & \ddots & \\
& \ddots & \ddots & \ddots \\
& & \ddots & T_{y,0} & T_{y,+1} \\
& & & T_{y,-1} & T_{y,0}
\end{pmatrix}_{n_{y}},
```

> "where the matrix elements depend on \beta_{x} through Eq. (6)."

**Figure 1** (the only figure in the fetched region; image file referenced as `2609.23523v1/F_1Dto2D.png`) — caption verbatim:

> "**Figure 1:** OBC energy spectra (black circles) for fixed L_{x}=100 and varying widths n_{y}=1,3,7, and 11. Red color intensity encodes the magnitude of the 1D GBZ deviation. Parameters: t_{x}=t_{-x}=i, t_{y}=t_{-y}=1/4. All computations use 16-digit numerical precision."

Exact numerical values for the Fig. 1 computation (verbatim from the caption, itemized):
- L_{x} = **100** (fixed)
- n_{y} ∈ {**1, 3, 7, 11**}
- t_{x} = t_{-x} = **i** (imaginary unit)
- t_{y} = t_{-y} = **1/4**
- all computations at **16-digit numerical precision**
- plot style: OBC energy spectra as **black circles** in the complex-energy plane; **red color intensity** encodes magnitude of the 1D GBZ deviation (i.e. |Δ(E)|, see Eq. 9).
- Derived (arithmetic, NOT quoted): total real-space OBC matrix sizes are L_x·n_y ∈ {100, 300, 700, 1100} → i.e. 100×100, 300×300, 700×700, 1100×1100. The paper does not state these sizes explicitly in the fetched region.

### 2.3 Section I.1 — the numerical protocol for "deviation from the 1D GBZ condition" (this is item (c) of the assignment)

> "Fig. 1 shows the OBC spectra for n_{y}=1,3,7, and 11 at fixed L_{x}=100. As n_{y} increases, the spectrum evolves from a curve toward an area-filling pattern, indicating the crossover from 1D to 2D. To quantify deviations from Eq. (4), we substitute each finite-size eigenvalue E into the ChE,"

**Eq. (8)** (numbered):
```
\det[E\mathbbm{I}_{n_{y}}-H_{n_{y}}(\beta_{x})]=0
```

> "and obtain **2n_{y} roots** \{\beta_{x,i}(E)\}, ordered by increasing modulus."

> "We define a quantity \Delta(E) to characterize the deviation from the 1D GBZ condition:"

**Eq. (9)** (numbered) — the deviation measure:
```
\Delta(E)=|\beta_{x,n_{y}p+1}(E)|-|\beta_{x,n_{y}p}(E)|.
```

> "Here **p=1** for this model. In Fig. 1, darker red indicates a larger deviation. As n_{y} increases, more OBC eigenvalues depart from the 1D GBZ condition, with a clear dependence on their position in the spectrum. The eigenvalues in the triangular regions near the imaginary axis nearly satisfy Eq. (4), whereas those in the left and right wings show larger deviations. These deviations increase with the transverse width."

Protocol as literally specified in the fetched text (step-by-step, with no additions):
1. Build the finite OBC Hamiltonian of the 2D model (Eq. 5) on an L_{x} × n_{y} rectangle with OBC in **both** directions; equivalently, the effective 1D chain of L_{x} supercells, each supercell carrying the n_{y}×n_{y} matrix Eq. (7) with β_{x}-dependent blocks Eq. (6).
2. Diagonalize to get finite-size eigenvalues E (spectrum plotted as black circles in Fig. 1).
3. For **each** finite-size eigenvalue E, substitute E into the ChE Eq. (8) and solve for the roots β_{x} — the polytope/root set \{\beta_{x,i}(E)\}, of cardinality **2n_{y}**.
4. Order those roots by increasing modulus.
5. Form Δ(E) = |β_{x, n_{y}p+1}(E)| − |β_{x, n_{y}p}(E)| with p = 1, i.e. Δ(E) = |β_{x, n_{y}+1}(E)| − |β_{x, n_{y}}(E)|.
6. Color each eigenvalue in the complex-energy plane by |Δ(E)| (darker red = larger deviation).

Consistency check (derived, not quoted): with the reduced model having n = n_y internal states and p = q = 1, the general count n_s = n(p+q) gives n_s = 2n_y — matching the "2n_{y} roots" statement, and the GBZ index in Eq. (4) is n_y·p = n_y (i.e. the n_y-th and (n_y+1)-th roots), matching Eq. (9).

**What the protocol does NOT specify in the fetched region (explicit negatives):**
- **No numerical method for the root solve is stated** — no companion matrix, no `numpy.roots`, no eigen-solver on a linearization, no polynomial-coefficient construction described, no library/toolkit named.
- **No tolerance / convergence criterion** is given for the finite-size eigenvalue computation, nor for judging "nearly satisfy Eq. (4)".
- **No definition of the color normalization** of Fig. 1 (linear/log/scale range, normalization constant) is given.
- **No tie-breaking rule** for equal-modulus roots (degenerate |β|) is given.
- **No explicit statement of the matrix sizes** used for the diagonalization (only L_x=100 and n_y values).
- **No pseudocode, no algorithm box, no code listing anywhere in Sections I–II.**
- The phrase "**16-digit numerical precision**" appears **only in the Fig. 1 caption**; it is not restated as a method description in the body, and it is not specified whether it refers to machine double precision or to a working-precision setting of the root finder.
- Number of figures in the fetched region: **1** (Fig. 1). The paper has 10 figures total (per arXiv comments), so 9 figures lie outside the fetched region.

### 2.4 Section I.1 — the three questions raised (item (d) of the assignment)

> "This example shows why the 1D problem must be examined before extending the theory to higher dimensions. It raises three questions:"
>
> "(i) For fixed n_{y}, does the spectrum converge to the 1D GBZ prediction as L_{x}\to\infty, or do deviations persist in this limit?"
>
> "(ii) If convergence occurs, how large must L_{x} be for the thermodynamic prediction to become reliable?"
>
> "(iii) If the required size is too large for direct numerical verification, how can the theory be tested numerically?"

> "These three questions will also be addressed in the final section."

**All questions raised in the fetched region, collected:**
- Q0 (Abstract / §I.1 motivating question): "how does this breakdown actually occur as a system gradually crosses over from 1D to 2D?" / "how does the 1D GBZ condition in Eq. (4) break down as the system crosses over from 1D to 2D, and how can this breakdown be characterized?"
- Q1 (§I, §2.1 above): "what determines \mu_{m,\rm GBZ}(k)?"
- Q2 (§I): "Understanding when and why the conventional condition fails is a central goal of this series."
- Q3 (§I.1): the three numbered questions (i)–(iii) above; per the text they are re-addressed "in the final section" (§VIII, outside the fetched region).
- Q4 (§I.2): the paper's own goal — "we review the 1D spectral approach, demonstrate that the existing theory cannot resolve these questions, and clarify where the fundamental difficulties lie."

### 2.5 Section I.2 — Organization of the paper (verbatim)

> "In this paper, we review the 1D spectral approach, demonstrate that the existing theory cannot resolve these questions, and clarify where the fundamental difficulties lie."

> "Section II formulates the 1D open boundary problem. Sections III and IV derive the boundary matrix equation and present exact finite-size calculations via **Vieta's formulas**. Sections V, VI, and VII review the spectral approach in the thermodynamic limit, its alternative formulations, and its limitations. Section VIII summarizes the results, discusses the open questions, and outlines **Papers II–IV**."

---

## 3. SECTION II — 1D Open Boundary Problem (readable portion only)

> "In the motivating example, we described a strip by grouping its n_{y} sites at each x coordinate into one supercell, reducing the problem to an effective 1D model. **We now consider the general 1D open boundary problem, which is also applicable to these quasi-1D systems.**"

> "Consider a general 1D n-band model whose non-Bloch Hamiltonian is a matrix Laurent polynomial of the complex variable \beta:"

**Eq. (10)** (numbered):
```
H(\beta)=\sum_{i=-p}^{q}T_{i}\beta^{i},
```

> "where p,q\in\mathbb{N}^{+} are the maximum hopping ranges to the right and left, respectively. The n\times n matrix T_{i} describes hopping from unit cell x+i to unit cell x. Here n counts the internal states per unit cell, including spin, orbital, sublattice, and layer degrees of freedom. The eigenenergy E of H(\beta) satisfies the ChE"

**Eq. (11)** (numbered):
```
f(E,\beta):=\det[E\mathbbm{I}_{n}-H(\beta)]=0.
```

> "When the longest-range hopping matrices satisfy the full-rank condition"

**Eq. (12)** (numbered):
```
\det[T_{-p}]\neq 0,\quad \det[T_{q}]\neq 0,
```

> "multiplying f(E,\beta) by \beta^{np} gives a polynomial of degree n_{s}\equiv n(p+q) with a nonzero constant term. Thus, for each E, the ChE has exactly n_{s} finite, nonzero roots, counted with multiplicity. **Unless stated otherwise, we assume the full-rank condition throughout.** An example with rank-deficient hopping matrices is the non-Hermitian Su–Schrieffer–Heeger (SSH) model, discussed in Sec. IV.2."

> "The hopping matrices \{T_{i}\} also determine the real-space Hamiltonian under OBCs. For a chain of L unit cells, with L\geq p+q so that the boundary regions do not overlap, the eigenvalue equation reads"

**Unnumbered block-matrix eigenvalue equation** — *this is exactly where the fetch truncates.* Reconstructable content from the readable LaTeX prefix (the matrix `\left[\begin{array}...\right]` acting on the column vector `[\psi(1);\psi(2);\vdots;\psi(p+1);\vdots;\psi(L-q);\vdots;\psi(L-1);\psi(L)]`, set equal to `E` times that same vector; the LaTeX contains `\hline\cr` separators, i.e. visible horizontal rules separating the first `p` cells and the last `q` cells from the bulk):

Block-Toeplitz structure as it appears:
- Row 1: `T_0  T_1  \cdots  T_q`
- Row 2: `T_{-1}  T_0  \ddots  T_{q-1}  T_q`
- Row 3: `\vdots  \ddots  \ddots  \ddots  \ddots  \ddots`
- (hline separator)
- Bulk rows: `T_{-p}  T_{-p+1}  \ddots  T_0  T_1  \ddots  T_q`; then `\ddots \ddots \ddots \ddots \ddots \ddots \ddots \ddots`; then `T_{-p} \ddots T_{-1} T_0 \ddots T_{q-1} T_q`
- (hline separator)
- Penultimate rows: `\ddots \ddots \ddots \ddots \ddots \vdots`; `T_{-p} T_{-p+1} \ddots T_0 T_1`
- Last row: `T_{-p} \cdots T_{-1} T_0`

Right-hand side vector (readable part): `[\psi(1); \psi(2); \vdots; \psi(p+1); \vdots; \psi(L-q); \vdots; \psi(L-1); \psi(L)]`, multiplied by `E`. The truncation occurs in the middle of the second copy of this vector inside the `\left[\begin{array}...\right]` (after `\psi(L)` there is a garbled `\hbox to0pt{...\lxSVG@begingroup@{\_sc` fragment, i.e. the alignment/hline machinery of the HTML renderer). **Everything after this point in Section II is not present in the fetched region.**

Notes on §II:
- The equation above carries **no visible equation number** in the fetched HTML.
- The condition `L\geq p+q` ("so that the boundary regions do not overlap") is the only stated finite-size constraint on L in the fetched region.
- No ansatz for ψ(x) (e.g. superposition of β^x solutions), no secular/boundary-matrix construction, no thermodynamic-limit statement is visible — all of that lies beyond the truncation.

---

## 4. Item (c) restated: exact numerical protocol present vs. absent

**Present in fetched region (verbatim-able):**
- Object computed: OBC energy spectra (eigenvalues of the finite real-space OBC Hamiltonian) for the model Eq. (5), Eq. (7); and, per eigenvalue, the roots of the ChE Eq. (8) (2n_y roots, ordered by increasing modulus).
- Deviation measure: Eq. (9), `\Delta(E)=|\beta_{x,n_{y}p+1}(E)|-|\beta_{x,n_{y}p}(E)|`, with `p=1`.
- The finite-size eigenvalue E used in Eq. (8) is explicitly a **finite-size** eigenvalue ("we substitute each finite-size eigenvalue E into the ChE").
- Parameters and precision: L_x=100, n_y=1,3,7,11, t_x=t_{-x}=i, t_y=t_{-y}=1/4, 16-digit numerical precision (caption of Fig. 1).
- Roots counted/ordered as in I.1 and II ("ordered by increasing modulus"; "counted with multiplicity").

**Absent from fetched region (must not be invented):**
- How the polynomial coefficients are assembled, how roots are obtained numerically (companion matrix / linearization / iterative), any library or version, any tolerance.
- Whether Δ(E) is used raw or normalized for the color map (the caption says only "Red color intensity encodes the magnitude of the 1D GBZ deviation", and the body says "darker red indicates a larger deviation").
- Any L_x-scaling study (the L_x→∞ question (i) is posed, not answered, in this region).
- Any statement of the actual root moduli values, any table of numbers, any convergence order, any error bars.
- Any statement of the total matrix dimension used for diagonalization.
- Any pseudocode/algorithm listing (explicitly: **no pseudocode in Sec. I–II**).
- No data-set/benchmark names (this is not an ML-style paper).

---

## 5. Item (e): code/data availability, acknowledgements, links (searched the whole fetched text)

- **Code availability statement: not present in fetched region.** No "Code available at", no GitHub/Zenodo/DOI-to-repository, no data availability sentence appears anywhere in the Abstract, §I, or the readable part of §II.
- **Data availability statement: not present in fetched region.**
- **Acknowledgements / funding statement: not present in fetched region.** (Acknowledgements are normally near the end of the paper, i.e. in §VIII, which is outside the fetched region.)
- **"Data availability" links: not present.** The only link-like artifacts in the fetched HTML are (a) arXiv's own boilerplate "Code, Data and Media Associated with this Article" widget block on the **abs page** (alphaXiv, CatalyzeX Code Finder, DagsHub, Gotit.pub, Hugging Face, ScienceCast, Replicate, TXYZ.AI) — these are arXiv platform widgets, **not** author-provided links, and none of them indicated an actual deposit for this paper; (b) author e-mail addresses `jphu@iphy.ac.cn` and `yangzs@xmu.edu.cn`; (c) normal hyperlinks to journal references and internal section anchors; (d) the license link "arXiv.org perpetual, non-exclusive license".
- **The only figure asset path visible** is the image file `2609.23523v1/F_1Dto2D.png` (Fig. 1).
- Access routes offered by arXiv for this paper: `/pdf/2609.23523v1`, `/html/2609.23523v1` (experimental), `/src/2609.23523v1` (TeX source). No companion-code link.

---

## 6. Numbered-equation inventory for the fetched region

| Eq. | Where | Content (cleaned) |
|---|---|---|
| (1) | §I | `\beta \equiv e^{ik} e^{\mu}` |
| (2) | §I | `\beta_{m,\rm GBZ}(k) \equiv e^{ik} e^{\mu_{m,\rm GBZ}(k)}` |
| (3) | §I | `\sigma_{\mathrm{OBC,bulk}}=\bigcup_{m}\{E_{m}(\beta_{m,\rm GBZ}(k)) \mid k\in[0,2\pi)\}` |
| (4) | §I | `|\beta_{np}(E)|=|\beta_{np+1}(E)|` (conventional 1D GBZ condition) |
| (5) | §I.1 | `H(\beta_x,\beta_y)=t_x\beta_x+t_{-x}/\beta_x+t_y\beta_x\beta_y+t_{-y}/(\beta_x\beta_y)` |
| (6) | §I.1 | `T_{y,0}=t_x\beta_x+t_{-x}/\beta_x`, `T_{y,+1}=t_y\beta_x`, `T_{y,-1}=t_{-y}/\beta_x` |
| (7) | §I.1 | `H_{n_y}(\beta_x)` = tridiagonal block-Toeplitz `n_y\times n_y` matrix with diagonal `T_{y,0}`, superdiagonal `T_{y,+1}`, subdiagonal `T_{y,-1}` |
| (8) | §I.1 | `\det[E\mathbbm{I}_{n_y}-H_{n_y}(\beta_x)]=0` → `2n_y` roots |
| (9) | §I.1 | `\Delta(E)=|\beta_{x,n_y p+1}(E)|-|\beta_{x,n_y p}(E)|`, `p=1` |
| — | §I | unnumbered: `|\beta_1(E)|\leq|\beta_2(E)|\leq\cdots\leq|\beta_{n_s}(E)|` |
| (10) | §II | `H(\beta)=\sum_{i=-p}^{q}T_i\beta^i` |
| (11) | §II | `f(E,\beta):=\det[E\mathbbm{I}_n-H(\beta)]=0` |
| (12) | §II | `\det[T_{-p}]\neq 0,\ \det[T_q]\neq 0` (full-rank condition) |
| — | §II | unnumbered block-matrix OBC eigenvalue equation (truncated mid-render) |

---

## 7. Titles of sections OUTSIDE the fetched region (from the HTML table of contents only — no body text read)

1. (Abstract) ✔ read
2. **I Introduction** ✔ read — I.1 Motivating example ✔, I.2 Organization of the paper ✔
3. **II 1D Open Boundary Problem** ✔ read up to the truncation point
4. III Boundary Matrix Equation — III.1 Boundary Matrix for a Finite Chain, III.2 Questions in the Thermodynamic Limit ✘
5. IV Exact Calculations for Finite Systems — IV.1 Example: Single-Band Model with Long-Range Hopping, IV.2 Example: Non-Hermitian SSH Model with Rank Deficiency ✘
6. V Spectral Approach I: Single-Band Case — V.1 GBZ Condition and Spectral Condensation, V.2 Predictions and Limitations of the GBZ ✘
7. VI Spectral Approach II: Alternative Formulations of the GBZ Condition — VI.0.1 **Root Criterion**, VI.0.2 **Amoeba Criterion**, VI.0.3 **Winding Number Criterion**, VI.0.4 **Spectral Potential Criterion**, VI.0.5 Summary and Scope of Validity ✘
8. VII Spectral Approach III: Limitations and Challenges — VII.1 Conventional and Anomalous Multi-Band GBZ Conditions, VII.2 Spectral Degeneracy Is Not Necessary, VII.3 Factorization of the ChE Is Not Sufficient ✘
9. VIII Summary and Discussion ✘
10. References ✘

(Reference hyperlink labels in the HTML run `bib.bib1` … `bib.bib48`, suggesting ~48 references; the bibliography text itself is outside the fetched region.)

---

## 8. Facts most relevant to a "brute-force SGBZ / amoeba" implementation (flagged for the parent, derived from the above, no invention)

- The 2D-to-effective-1D reduction used in §I.1 yields a **single-band 2D model whose reduced problem is an n_y-band 1D problem** (n = n_y, p = q = 1), so the ChE in β_x is a polynomial with **n_s = n(p+q) = 2n_y** finite nonzero roots and the GBZ index pair is (n_y, n_y+1) — exactly the pair Eq. (9) differences. This gives a concrete, small, exactly controllable testbed: H(β) = t_xβ + t_{-x}/β + t_yββ_y + t_{-y}/(ββ_y) with t_x=t_{-x}=i, t_y=t_{-y}=1/4.
- The paper's own definition of the "deviation" quantity is **one-sided and signed** (a difference of moduli, not an absolute value inside the definition; absolute value only in the color-magnitude description), and ordering is **by increasing modulus with multiplicity counted** — relevant to how a brute-force root/ordering routine must break ties.
- The paper explicitly asserts the full-rank condition `det[T_{-p}]≠0, det[T_q]≠0` as the standing assumption and defers the rank-deficient (non-Hermitian SSH) case to §IV.2.
- §VI.0.2 "Amoeba Criterion" is the paper's amoeba-based formulation of the GBZ condition — outside my assigned/readable region, but it exists and is titled exactly that.

---

## 9. Integrity statement

Everything above is transcribed from the single fetched HTML region (Abstract, §I, §I.1, §I.2, §II up to the truncation) plus arXiv metadata from https://arxiv.org/abs/2609.23523v1. Nothing was inferred about unread text; every absence is marked "not present in fetched region". The truncation is real and reproducible: three separate fetch attempts (plain URL, `#S2` fragment, and after the abs page) all cut at the identical point inside the unnumbered block-matrix equation of §II, and the PDF route is unavailable to this tool.

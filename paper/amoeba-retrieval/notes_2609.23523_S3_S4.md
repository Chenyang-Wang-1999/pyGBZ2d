# Extraction attempt: arXiv:2609.23523v1 — Sections III and IV

**Paper:** Zeqi Xu, Jiangping Hu, Zhesen Yang, "Universal Generalized Brillouin Zone Theory I: Review of the Spectral Approach", arXiv:2609.23523v1 [cond-mat.mes-hall], submitted Sun, 20 Sep 2026 10:18:51 UTC, comments field "16 pages, 10 figures", DOI 10.48550/arXiv.2609.23523.

**Assigned reading region = Section III ("Boundary Matrix Equation", incl. III.1 "Boundary Matrix for a Finite Chain", III.2 "Questions in the Thermodynamic Limit") + Section IV ("Exact Calculations for Finite Systems", incl. IV.1 "Example: Single-Band Model with Long-Range Hopping", IV.2 "Example: Non-Hermitian SSH Model with Rank Deficiency").**

## 0. HEADLINE VERDICT — READ THIS FIRST

**Sections III and IV were NOT retrieved. Not one line of their body text, not one equation (numbered ≥ 13), not one figure caption, and not one algorithm/pseudocode/numerics statement from them is in my possession.**

This is a hard tool-level limitation, not a partial-read. The `web_fetch` tool available in this session truncates the arXiv HTML page at a fixed character budget, and that budget is exhausted **inside Section II**, before Section III begins. Every retrieval route that reaches the full text is blocked (details in §1). A background subagent was dispatched to sweep ~15 alternative endpoints and independently confirmed the identical truncation point, with zero III/IV body text obtained.

Therefore: **every item in the requested extraction list (1)–(7) is marked "not retrieved" below**, with an explicit statement of what I could and could not see. Nothing has been inferred, reconstructed, or paraphrased into the gaps.

Two sibling workspace files record the same limitation from other angles and should be read alongside this one:
- `notes_2609.23523_S1_S2.md` — sibling agent's extraction for Abstract/I/II; its §7 marks Sections III–VIII as unread (lines 322–323).
- `notes_2609.23523_endpoint-retrieval-report.md` — the background subagent's per-URL retrieval log.

---

## 1. RETRIEVAL ATTEMPT LOG (what worked, what failed, and exactly where it truncates)

### 1.1 The truncation is real, reproducible, and identical everywhere

`https://arxiv.org/html/2609.23523v1` returns HTTP 200 and truncates at a **fixed point inside Section II**, immediately after the sentence:

> "The hopping matrices $\{T_i\}$ also determine the real-space Hamiltonian under OBCs. For a chain of $L$ unit cells, with $L \geq p+q$ so that the boundary regions do not overlap, the eigenvalue equation reads"

— i.e. at the first large $nL \times nL$ block-matrix display (the unnumbered equation that follows Eq. (12) in the paper's numbering). The tool emits: *"(Content truncated. Fetch a more specific URL or section for the full text.)"*

The **last characters actually received**, verbatim (identical on every successful endpoint):

```
\end{array}\right\]\left\[\begin{array}\[\]{c}\psi(1)\\\\[6.45831pt\] \psi(2)\\ \vdots\\ \hline\cr\\\[-6.45831pt\] \psi(p{+}1)\\\\[2.15277pt\] \vdots\\\[4.30554pt\] \psi(L{-}q)\\\\[2.15277pt\] \hline\cr\\\[-10.76385pt\] \vdots\\\[5.16663pt\] \psi(L{-}1)\\ \\\[-5.16663pt\] \psi(L)\\ \end{array}\right\]=E\left\[\begin{array}\[\]{c}\hbox to0pt{\vbox to0pt{\pgfpicture\makeatletter\hbox{\hskip 0.0pt\lower 0.0pt\hbox to0.0pt{\lxSVG@begingroup@{\_sc
```

So the readable region is: **HTML navigation/TOC (all ten section headings, titles only) + title/authors/abstract + Section I (all, including I.1 and I.2) + Section II up to and including Eqs. (10), (11), (12) and the sentence introducing the block-matrix display.** Nothing after that point.

### 1.2 Per-URL results (all attempted in this session)

| URL | Result |
|---|---|
| `https://arxiv.org/html/2609.23523v1` | HTTP 200, truncated at the §II block-matrix display |
| `https://export.arxiv.org/html/2609.23523v1` | HTTP 200, **byte-identical** content and truncation point |
| `https://arxiv.org/html/2609.23523` | HTTP 200, byte-identical (redirects to v1) |
| `https://arxiv.org/html/2609.23523v1/` (trailing slash) | HTTP 200, byte-identical |
| `https://arxiv.org/html/2609.23523v1#S3` | HTTP 200, byte-identical — fragments are **not** honoured by the fetcher |
| `.../#S4`, `.../#S4.SS1` | HTTP 200, byte-identical — fragments ignored |
| `https://arxiv.org/html/2609.23523v1?output=1` | HTTP 200, byte-identical — no-op |
| `https://arxiv.org/html/2609.23523v1?full=1` | HTTP 200, byte-identical — no-op |
| `https://ar5iv.labs.arxiv.org/html/2609.23523` | **Error:** "cross-origin redirect to https://arxiv.org is not followed automatically" (retry against arxiv.org gives the abstract page only) |
| `https://ar5iv.org/abs/2609.23523` | **Error:** cross-origin redirect chain, dead end |
| `https://www.arxiv-vanity.com/papers/2609.23523/` | **Error:** cross-origin redirect to ar5iv, dead end |
| `https://arxiv.org/pdf/2609.23523v1` | **Error:** `unsupported content type "application/pdf"` — this tool cannot read PDFs |
| `https://arxiv.org/pdf/2609.23523v1?download=true` | same PDF content-type refusal |
| `https://arxiv.org/e-print/2609.23523v1` (TeX source) | **Error:** `unsupported content type "application/gzip"` (per subagent) — cannot untar |
| `https://arxiv.org/html/2609.23523v1/2609.23523v1.pdf` | HTTP 404 |
| `https://arxiv.org/html/2609.23523v1/index.html` | HTTP 404 ("No HTML for '2609.23523v1'") |
| `https://arxiv.org/html/2609.23523v2` | HTTP 404 — **no v2 exists**; v1 is the only version |
| `https://www.arxiv.org/html/2609.23523v1` | **Error:** cross-origin redirect |
| `https://web3.arxiv.org/html/2609.23523v1` | **Error:** `getaddrinfo ENOTFOUND` |
| `https://r.jina.ai/https://arxiv.org/html/2609.23523v1` | **Error:** `TypeError: fetch failed` |
| `https://huggingface.co/papers/2609.23523` | **Error:** `TypeError: fetch failed` |
| `https://www.alphaxiv.org/abs/2609.23523` | HTTP 200 — **abstract + metadata only**, no body |
| `https://www.alphaxiv.org/pdf/2609.23523` | HTTP 200 — client-rendered shell, no body text |
| `https://www.semanticscholar.org/arxiv/2609.23523` | HTTP 404 |
| `https://api.semanticscholar.org/graph/v1/paper/arXiv:2609.23523` | HTTP 429 (rate limited) |
| `https://api.openalex.org/works/doi:10.48550/arXiv.2609.23523` | **Error:** `TypeError: fetch failed` |
| `https://openalex.org/works?search=2609.23523` | HTTP 403 (Cloudflare check) |
| `https://arxiv-org.ezproxy.obspm.fr/abs/2609.23523` | **Error:** cross-origin redirect to login |
| `https://arxiv-org.ezproxy.obspm.fr/html/2609.23523v1` | **Error:** (http and https) `TypeError: fetch failed` |
| `https://search.marcia.arxiv.org/abs/2609.23523` | **Error:** `getaddrinfo ENOTFOUND` |
| `https://arxiv.org/a/xu_z_1.html` | HTTP 404 |
| `https://arxiv.org/abs/2609.23523v1` | HTTP 200 — **abstract page only** (metadata, no body) |
| `https://export.arxiv.org/api/query?id_list=2609.23523` | HTTP 200 — Atom metadata only (title, authors, abstract, comment "16 pages, 10 figures"); **no body text** |

Additional mirrors tried by the background subagent (all without body text): `emergentmind.com/papers/2609.23523` (abstract only), `papers.cool/arxiv/2609.23523` (abstract + metadata), `hjfy.top/arxiv/2609.23523` (JS shell), `themoonlight.io/en/review/...` (HTTP 429 bot checkpoint), `textise.net`, `synthical.com` (fetch failed). arXiv HTML is a single monolithic LaTeXML page, so there is no section-split endpoint, and the inline-SVG math fallback of the §II matrix consumes the remaining budget.

### 1.3 What the truncation means, honestly stated

- I saw **Sections I and II only** (plus the TOC headings of III–VIII and the References heading).
- **I did not see any of Section III or Section IV.** I have no body text, no equations, no figure captions, no numerical recipe, no tolerances, no code statement from them.
- I therefore **cannot claim** to know how the boundary matrix equation is defined, how Vieta's formulas are used, what polynomial degrees arise in IV.1/IV.2, what $L$/$n$/$p$/$q$ values were used, what matrix sizes were diagonalized, what numerical precision or tolerances were chosen, what failure modes were reported, or whether code/data is available.
- **All of those fields below are recorded as "not retrieved".** No field has been filled with an inference dressed as a reading.

---

## 2. THE ASSIGNED ITEMS, FIELD BY FIELD

### Item 1 — The boundary matrix equation (definition, size, conversion of the finite-size OBC problem into an algebraic condition)

**NOT RETRIEVED.** I did not see a single line of Section III.

From the TOC only, I can report the section's existence and internal structure (`#S3`, `#S3.SS1`, `#S3.SS2`):

> "III Boundary Matrix Equation — III.1 Boundary Matrix for a Finite Chain; III.2 Questions in the Thermodynamic Limit"

Cross-reference from Section I.2 (the only description of Section III anywhere in the readable region):

> "Sections III and IV derive the boundary matrix equation and present exact finite-size calculations via Vieta's formulas."

That single sentence is the **entirety** of what I can honestly report about Section III's content. It confirms only that (a) a "boundary matrix equation" exists, (b) it is *derived* in Sections III–IV jointly, and (c) Vieta's formulas are involved in the finite-size calculation. It gives no equation, no size, no degree, no construction.

Related machinery that IS readable, but belongs to Section II (so it is the *precondition*, not the boundary matrix equation itself):

- Eq. (10), §II: `H(\beta)=\sum_{i=-p}^{q}T_{i}\beta^{i}` — the general 1D $n$-band non-Bloch Hamiltonian as a matrix Laurent polynomial; $p,q\in\mathbb{N}^{+}$ "the maximum hopping ranges to the right and left, respectively"; $T_i$ is $n\times n$ and "describes hopping from unit cell $x+i$ to unit cell $x$".
- Eq. (11), §II: `f(E,\beta):=\det[E\mathbbm{I}_{n}-H(\beta)]=0` — the characteristic equation (ChE).
- Eq. (12), §II (full-rank condition): `\det[T_{-p}]\neq 0,\quad \det[T_{q}]\neq 0`.
- Verbatim: *"multiplying $f(E,\beta)$ by $\beta^{np}$ gives a polynomial of degree $n_{s}\equiv n(p+q)$ with a nonzero constant term. Thus, for each $E$, the ChE has exactly $n_{s}$ finite, nonzero roots, counted with multiplicity. Unless stated otherwise, we assume the full-rank condition throughout. An example with rank-deficient hopping matrices is the non-Hermitian Su–Schrieffer–Heeger (SSH) model, discussed in Sec. IV.2."*
- Verbatim, the sentence immediately before the truncation: *"The hopping matrices $\{T_{i}\}$ also determine the real-space Hamiltonian under OBCs. For a chain of $L$ unit cells, with $L\geq p+q$ so that the boundary regions do not overlap, the eigenvalue equation reads"* → then the $nL\times nL$ block-Toeplitz matrix of blocks $T_{-p},\dots,T_q$, multiplied into the column vector $[\psi(1);\psi(2);\dots;\psi(p+1);\dots;\psi(L-q);\dots;\psi(L-1);\psi(L)]$ and set equal to $E$ times that same vector. The LaTeX contains `\hline\cr` separators marking the first $p$ cells and the last $q$ cells (the boundary regions) off from the bulk. **The render truncates mid-vector and I never see anything after it.**

I record explicitly that I have **no** information on: the size of the boundary matrix, whether it is $np\times np$, $nq\times nq$, $(np+nq)\times(np+nq)$ or otherwise; how the boundary conditions are imposed (secular determinant? rank condition on a transfer/companion matrix? vanishing of a submatrix?); or how that condition is claimed to become polynomial in $E$.

### Item 2 — Deterministic finite-size spectra via Vieta's formulas (companion matrix / resultant / characteristic polynomial in $E$ and $\beta$; polynomial degrees; closed form or root-finding)

**NOT RETRIEVED.**

The only evidence in the readable region that such a procedure exists is the §I.2 sentence: *"Sections III and IV derive the boundary matrix equation and present exact finite-size calculations via Vieta's formulas."*

Not retrieved, explicitly:
- whether a companion matrix is formed;
- whether a resultant (e.g. $\mathrm{Res}_\beta(f(E,\beta), g(E,\beta))=0$) is computed;
- whether a characteristic polynomial in $E$ alone is assembled, and of what degree;
- the degree and the $E$-dependence of any polynomial in $\beta$;
- any closed-form expression for finite-size $E$;
- any root-finding procedure (Newton, Aberth, Ehrlich–Aberth, Weierstrass/Durand–Kerner, eigenvalue-of-companion, etc.).

Guessing any of these would be fabrication, so I do not.

### Item 3 — Algorithm steps, pseudocode, or explicit numerical recipe

**NOT RETRIEVED for Sections III–IV.** No pseudocode box, algorithm listing, or numbered recipe is present anywhere in the readable region (Abstract, §I, §II) either — I searched the received text for it. The paper's own algorithmic content for III/IV is therefore entirely unread.

The nearest thing to a *numerical recipe* that I can honestly report is the **Section I.1 deviation-measurement protocol**, which is readable and is the only explicit step-by-step procedure in the text I actually saw. It is quoted here because it is the paper's one fully visible numerical pipeline, and it is the protocol that any brute-force GBZ/amoeba implementation in this workspace would be reproducing. Verbatim anchors:

> "To quantify deviations from Eq. (4), we substitute each finite-size eigenvalue $E$ into the ChE,"

**Eq. (8)** (§I.1): `\det[E\mathbbm{I}_{n_{y}}-H_{n_{y}}(\beta_{x})]=0`

> "and obtain $2n_{y}$ roots $\{\beta_{x,i}(E)\}$, ordered by increasing modulus."

**Eq. (9)** (§I.1): `\Delta(E)=|\beta_{x,n_{y}p+1}(E)|-|\beta_{x,n_{y}p}(E)|`

> "Here $p=1$ for this model. In Fig. 1, darker red indicates a larger deviation. As $n_{y}$ increases, more OBC eigenvalues depart from the 1D GBZ condition, with a clear dependence on their position in the spectrum. The eigenvalues in the triangular regions near the imaginary axis nearly satisfy Eq. (4), whereas those in the left and right wings show larger deviations. These deviations increase with the transverse width."

Restated as the steps the text literally specifies (nothing added):
1. Diagonalize the finite OBC Hamiltonian of the 2D model (Eq. 5) on an $L_x\times n_y$ rectangle with OBC in **both** directions; equivalently an effective 1D chain of $L_x$ supercells, each carrying the $n_y\times n_y$ matrix Eq. (7) with $\beta_x$-dependent blocks Eq. (6).
2. Take each finite-size eigenvalue $E$ from that diagonalization.
3. Substitute $E$ into the ChE Eq. (8) and solve for $\beta_x$ → a set $\{\beta_{x,i}(E)\}$ of cardinality $2n_y$.
4. Order those roots by increasing modulus ("counted with multiplicity" per §II).
5. Form $\Delta(E)=|\beta_{x,n_y+1}(E)|-|\beta_{x,n_y}(E)|$ (since $p=1$).
6. Colour each eigenvalue in the complex-energy plane by the magnitude of the deviation.

Consistency check (derived, not quoted): the reduced model has $n=n_y$ internal states and $p=q=1$, so $n_s=n(p+q)=2n_y$ — matching the stated "2n_y roots", and the conventional GBZ pair $(np, np+1)=(n_y,n_y+1)$ matches Eq. (9). This corroborates the §II counting rule; it does not supply the missing Section III/IV algorithm.

### Item 4 — Accuracy/precision control (digits, root ordering by modulus, degeneracy/clustering handling, tolerances, verification against direct diagonalization)

**NOT RETRIEVED for Sections III–IV.** Explicit negatives for the whole readable region:

- The only precision statement anywhere in the text I received is inside the Fig. 1 caption (§I.1): *"All computations use 16-digit numerical precision."* It is **not** restated in the body, and it is **not** specified whether it means IEEE-754 double (machine) precision, or a working-precision setting of a root finder (e.g. `mpmath.mp.dps = 16`). It also says nothing about the III/IV finite-size calculations.
- **Root ordering:** stated in §I.1 and §II as "ordered by increasing modulus" / "counted with multiplicity". No tie-breaking rule for equal-modulus roots is given, and no criterion for how to order degenerate $|\beta|$ is given.
- **Degeneracy / root-clustering handling: not stated** anywhere in the readable region.
- **Tolerances / convergence criteria: not stated.**
- **Verification against direct diagonalization:** the §I.1 protocol *does* use direct diagonalization of the finite OBC matrix to generate the $E$ values that are then fed into the ChE, so the ChE root solve is exercised *at* the finite-size eigenvalues; but no quantitative cross-check, residual tolerance, or error estimate is reported.
- No error bars, no convergence order, no residual norms, no tables of numbers appear in the readable region.

### Item 5 — Worked examples (exact models, hopping parameters, chain lengths $L$, internal states $n$, ranges $p$ and $q$, matrix/polynomial sizes, figure numbers and captions)

**NOT RETRIEVED for Section IV.1 and Section IV.2**, which are precisely the two requested worked examples. I saw no part of either.

- **IV.1 "Example: Single-Band Model with Long-Range Hopping"** — TOC heading only (`#S4.SS1`). No model, no parameters, no figure, no caption, no $L$, no $p,q$ values, no polynomial degree. **Not retrieved.**
- **IV.2 "Example: Non-Hermitian SSH Model with Rank Deficiency"** — TOC heading only (`#S4.SS2`), plus one pointer to it from §II: *"An example with rank-deficient hopping matrices is the non-Hermitian Su–Schrieffer–Heeger (SSH) model, discussed in Sec. IV.2."* That pointer establishes only that the non-Hermitian SSH model is the rank-deficient example (i.e. $\det[T_{-p}]=0$ and/or $\det[T_q]=0$ violates Eq. (12)). **No model details, no parameters, no figures, no captions retrieved.**

The **only** worked numerical example visible anywhere in my readable region is **Fig. 1 in §I.1** (the 2D→1D crossover motivating example). Its caption, verbatim:

> "**Figure 1:** OBC energy spectra (black circles) for fixed $L_{x}=100$ and varying widths $n_{y}=1,3,7$, and $11$. Red color intensity encodes the magnitude of the 1D GBZ deviation. Parameters: $t_{x}=t_{-x}=i$, $t_{y}=t_{-y}=1/4$. All computations use 16-digit numerical precision."

Itemised values from that caption (verbatim content, not inference):
- $L_x = 100$ (fixed); $n_y \in \{1, 3, 7, 11\}$
- $t_x = t_{-x} = i$ (imaginary unit); $t_y = t_{-y} = 1/4$
- 16-digit numerical precision
- black circles = OBC energy spectra in the complex-energy plane; red colour intensity = magnitude of the 1D GBZ deviation

Associated model, readable in §I.1:

**Eq. (5):** `H(\beta_{x},\beta_{y})=t_{x}\beta_{x}+t_{-x}/\beta_{x}+t_{y}\beta_{x}\beta_{y}+t_{-y}/(\beta_{x}\beta_{y})`

**Eq. (6):** `T_{y,0}(\beta_{x})=t_{x}\beta_{x}+t_{-x}/\beta_{x}`, `T_{y,+1}(\beta_{x})=t_{y}\beta_{x}`, `T_{y,-1}(\beta_{x})=t_{-y}/\beta_{x}`

**Eq. (7):** `H_{n_{y}}(\beta_{x})` = the $n_y\times n_y$ tridiagonal block-Toeplitz matrix with diagonal $T_{y,0}$, superdiagonal $T_{y,+1}$, subdiagonal $T_{y,-1}$ (open boundaries in $y$).

Derived (arithmetic, **not** stated in the paper): reduced model has $n=n_y$, $p=q=1$, so the total real-space OBC matrix is $L_x n_y\times L_x n_y$, i.e. $100\times100$, $300\times300$, $700\times700$, $1100\times1100$ for $n_y=1,3,7,11$; and the $\beta_x$-polynomial in Eq. (8) has $n_s=n(p+q)=2n_y$ finite nonzero roots (2, 6, 14, 22). The paper does not print these sizes explicitly in the readable region.

**Figure count:** 1 figure (Fig. 1) lies in the readable region. Per the arXiv comments field the paper has **10 figures**, so **9 figures lie outside what I could read** — including every figure in Sections III, IV, V, VI, VII, VIII. Their numbers and captions are **not retrieved**.

### Item 6 — Known failure modes of the exact finite-size calculations (numerical instability, ill-conditioning of the polynomial, etc.)

**NOT RETRIEVED.** No statement about conditioning, instability, overflow/underflow, root clustering, spurious roots, or loss of precision in the finite-size calculation appears anywhere in my readable region.

What the readable region *does* contain, and what a reviewer might confuse for this item, is a different kind of admitted limitation — the paper's three open questions about the *1D GBZ theory itself*, raised at the end of §I.1 and promised to be revisited in §VIII (outside the readable region):

> "(i) For fixed $n_{y}$, does the spectrum converge to the 1D GBZ prediction as $L_{x}\to\infty$, or do deviations persist in this limit?"
> "(ii) If convergence occurs, how large must $L_{x}$ be for the thermodynamic prediction to become reliable?"
> "(iii) If the required size is too large for direct numerical verification, how can the theory be tested numerically?"

These are questions about theory validity and required system size, **not** numerical failure modes of the exact finite-size method. I record them as such and do not present them as Item 6 content.

### Item 7 — Code/data availability statement

**NOT RETRIEVED** — and specifically, **no code or data availability statement appears in the Abstract, §I, §I.1, §I.2, or the readable part of §II.** Statements of this kind normally sit in the acknowledgements/final section (§VIII), which is outside the readable region, so their absence here is not evidence that none exists.

What is verifiable: the arXiv **abs page** for this paper carries only arXiv's own generic "Code, Data and Media Associated with this Article" platform widgets (alphaXiv, CatalyzeX Code Finder, DagsHub, Gotit.pub, Hugging Face, ScienceCast, Replicate, TXYZ.AI) — these are arXiv platform boilerplate, not author deposits, and none of them indicated an actual deposit for this paper. Access routes arXiv itself lists are only: PDF (`/pdf/2609.23523v1`), HTML experimental (`/html/2609.23523v1`), and TeX source (`/src/2609.23523v1`). No companion-repository link was found.

---

## 3. WHAT IS VERBATIM-AVAILABLE FOR SECTION III/IV CONTEXT (all of it)

This is the complete set of sentences anywhere in my readable region that refer to Sections III or IV. Nothing else exists.

1. From §I.2 (Organization of the paper):
   > "Section II formulates the 1D open boundary problem. Sections III and IV derive the boundary matrix equation and present exact finite-size calculations via Vieta's formulas. Sections V, VI, and VII review the spectral approach in the thermodynamic limit, its alternative formulations, and its limitations. Section VIII summarizes the results, discusses the open questions, and outlines Papers II–IV."
2. From §II:
   > "An example with rank-deficient hopping matrices is the non-Hermitian Su–Schrieffer–Heeger (SSH) model, discussed in Sec. IV.2."
3. HTML table of contents / navigation (headings and anchor ids only, no body):
   - `#S3` — "III Boundary Matrix Equation"
   - `#S3.SS1` — "III.1 Boundary Matrix for a Finite Chain"
   - `#S3.SS2` — "III.2 Questions in the Thermodynamic Limit"
   - `#S4` — "IV Exact Calculations for Finite Systems"
   - `#S4.SS1` — "IV.1 Example: Single-Band Model with Long-Range Hopping"
   - `#S4.SS2` — "IV.2 Example: Non-Hermitian SSH Model with Rank Deficiency"

Also worth flagging for the workspace, as a *pointer* rather than content: §VI.0.2 of the TOC is titled **"Amoeba Criterion"** (`#S6.SS0.SSS2`), and §VI.0.1/VI.0.3/VI.0.4 are "Root Criterion", "Winding Number Criterion", "Spectral Potential Criterion". These sections define this paper's alternative formulations of the GBZ condition — directly relevant to the brute-force SGBZ/amoeba work in this repository — but their body text is outside the readable region and was **not retrieved**.

---

## 4. NUMBERED-EQUATION INVENTORY OF WHAT I ACTUALLY SAW

| Eq. | Section | Content (cleaned LaTeX) |
|---|---|---|
| (1) | §I | `\beta \equiv e^{ik}e^{\mu}` |
| (2) | §I | `\beta_{m,\rm GBZ}(k)\equiv e^{ik}e^{\mu_{m,\rm GBZ}(k)}` |
| (3) | §I | `\sigma_{\mathrm{OBC,bulk}}=\bigcup_{m}\{E_{m}(\beta_{m,\rm GBZ}(k))\mid k\in[0,2\pi)\}` |
| (4) | §I | `|\beta_{np}(E)|=|\beta_{np+1}(E)|` (conventional 1D GBZ condition) |
| (5) | §I.1 | `H(\beta_x,\beta_y)=t_x\beta_x+t_{-x}/\beta_x+t_y\beta_x\beta_y+t_{-y}/(\beta_x\beta_y)` |
| (6) | §I.1 | `T_{y,0}=t_x\beta_x+t_{-x}/\beta_x`; `T_{y,+1}=t_y\beta_x`; `T_{y,-1}=t_{-y}/\beta_x` |
| (7) | §I.1 | `H_{n_y}(\beta_x)`: $n_y\times n_y$ tridiagonal block-Toeplitz, diag $T_{y,0}$, superdiag $T_{y,+1}$, subdiag $T_{y,-1}$ |
| (8) | §I.1 | `\det[E\mathbbm{I}_{n_y}-H_{n_y}(\beta_x)]=0` → $2n_y$ roots |
| (9) | §I.1 | `\Delta(E)=|\beta_{x,n_y p+1}(E)|-|\beta_{x,n_y p}(E)|`, $p=1$ |
| — | §I | unnumbered: `|\beta_1(E)|\leq|\beta_2(E)|\leq\cdots\leq|\beta_{n_s}(E)|` |
| (10) | §II | `H(\beta)=\sum_{i=-p}^{q}T_i\beta^i` |
| (11) | §II | `f(E,\beta):=\det[E\mathbbm{I}_n-H(\beta)]=0` |
| (12) | §II | `\det[T_{-p}]\neq 0,\ \det[T_q]\neq 0` (full-rank condition) |
| — | §II | unnumbered $nL\times nL$ block-Toeplitz OBC eigenvalue equation — **truncation occurs mid-render here** |

**Equations (13) onward were never seen.** Since Section III is where the boundary matrix equation is derived, the paper's numbering for it is unknown to me.

---

## 5. HOW TO UNBLOCK THIS (concrete, for the parent agent)

The content is not unobtainable — it is unobtainable **through this session's `web_fetch`**. Any of the following would work, and I can read local files normally:

1. **Fetch the PDF or TeX source outside the fetch-tool path** (e.g. `curl -L -o p.pdf https://arxiv.org/pdf/2609.23523v1` or `curl -L -o src.tar.gz https://arxiv.org/e-print/2609.23523v1`, then `pdftotext p.pdf` / `tar xzf src.tar.gz`), place the extracted text in this workspace, and hand me the path. I will then extract Items 1–7 properly.
2. **Paste the Section III and IV text** (from the PDF or the TeX source) into a file in the workspace and point me at it.
3. Find a mirror that **splits sections** or **strips the inline-SVG math fallback** (the multi-kB `\lxSVG@...` blob in the §II matrix display is what exhausts the budget — a mirror without it would likely fit far more).
4. A fetch tool supporting **byte offsets or per-section extraction** would also solve it.

Until one of those happens, Sections III–IV of this paper should be treated in the literature-survey pipeline as **not read**, not as read-and-summarized.

---

## 6. INTEGRITY STATEMENT

Every statement above is either (a) a verbatim quote from the received text, (b) a clean-up of the doubled MathJax escaping inside a quoted equation (a pure notational fix that adds no content), (c) arXiv metadata read directly from `https://arxiv.org/abs/2609.23523v1` and the arXiv API, or (d) explicitly labelled as derived arithmetic. The retrieval failures in §1 are reported as observed. **No content of Sections III or IV — no equation, no parameter, no figure, no algorithm step, no tolerance, no failure mode, no code statement — has been invented, inferred, or reconstructed to fill the gap.** Every requested field that could not be read is marked "NOT RETRIEVED" in capital letters. Independent confirmation of the truncation boundary came from a separately-dispatched subagent sweep over ~15 additional endpoints, logged in `notes_2609.23523_endpoint-retrieval-report.md`.

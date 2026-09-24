# Derivative papers: how later work states and quotes Wang–Song–Wang (WSW)

**Target work (WSW):** Hong-Yi Wang, Fei Song, Zhong Wang, "Amoeba Formulation of Non-Bloch Band Theory in Arbitrary Dimensions", arXiv:2212.11743v3, Phys. Rev. X **14**, 021011 (2024).

**Purpose of this file:** verbatim extracts from *derivative* papers (papers that cite/restate/extend WSW), each labelled with the derivative paper it comes from and the exact URL it was retrieved from. Nothing here is quoted from WSW itself unless explicitly marked as such (see Appendix E). All text below was copied from tool output in this session; anything I could not read is listed at the end.

**Retrieval routes — what worked and what did not (verified this session):**
- `arxiv.org/html/<id>` → usable, but **truncated at roughly the end of Section II** for long papers, because the LaTeXML render duplicates every display equation (LaTeX source + MathML), roughly tripling the character count.
- `ar5iv.labs.arxiv.org/html/<id>` → **reaches noticeably further** than `arxiv.org/html` for the same paper (it also exposes author–year inline citations, which lets reference numbers be mapped to names). Best general route found. Note it may render an *older version* than `arxiv.org/html`.
- `arxiv.org/pdf/...`, `harvest.aps.org/.../fulltext`, `nature.com/...pdf` → rejected by the fetch tool: `unsupported content type "application/pdf"`. **No PDF route is available.**
- `nature.com/articles/s42005-026-02546-2` → blocked (`cross-origin redirect to https://idp.nature.com`); `link.springer.com` mirror → `fetch failed`; Europe PMC REST lookup by DOI → `"hitCount":0`.
- `web_search` → its index **does** contain the full text of these papers, and returns verbatim in-document sentences as "sources" with fragment anchors (`#page#para` for PDFs). Used for the fragments that lie beyond the fetch truncation point. Snippets are usually one sentence, truncated at ~150 characters.

---

## A. Kaneshiro & Peters — two papers by the same authors

### A.1 arXiv:2511.11349 — "Wiener-Hopf factorization and non-Hermitian topology for Amoeba formulation in one-dimensional multiband systems" (Phys. Rev. Research **8**, 013292)

Sources used: `https://arxiv.org/html/2511.11349v2`, `https://arxiv.org/html/2511.11349v1`, `https://ar5iv.labs.arxiv.org/html/2511.11349`.

**(a) The OBC-spectrum criterion in amoeba/Ronkin terms (Section II.2, "Amoeba formulation for class A").** After deriving the piecewise-linear Ronkin function, the paper states (verbatim, v2 HTML, Eq. (17) region):

> "Using this compact expression, one finds that Rσ(μ) attains its minimum within the interval [μp, μp+1], which precisely corresponds to the GBZ condition [7],[9]. In this region, the transformed symbol σμ is topologically trivial, corresponding to a vanishing winding number W[σμ]=0. For energies on the OBC spectrum, this interval collapses to a single value of μ, whose magnitude gives the inverse localization length of the corresponding eigenstate."

The same section defines the objects it uses (verbatim, Eqs. (14)–(15), v2 HTML):

> "This discussion is equivalent to the optimization problem
>   φ(E) = min_μ Rσ(μ),   (14)
> where Rσ is the Ronkin function defined as
>   Rσ(μ) = ∫₀^{2π} dk/2π ln|det σ_μ(e^{ik})|,   (15)
> and its derivative with respect to μ directly yields the winding number, W[σμ]."

and (verbatim):

> "This theorem states that the OBC and PBC spectral potentials coincide for such E." (i.e. when the winding number W[σ] in their Eq. (10) is "well-defined and vanishes")

*Note:* the three formulations in the task brief (no central hole / minimum attained on a plateau of the Ronkin function / vanishing winding number) appear in this paper as: (i) the optimum sitting in the interval [μ_p, μ_{p+1}] where "the transformed symbol σ_μ is topologically trivial, corresponding to a vanishing winding number W[σμ] = 0"; (ii) "For energies on the OBC spectrum, this interval collapses to a single value of μ" (the plateau-collapse form). The explicit words "central hole" / "hole closing" are used in the companion paper A.2, not here.

**(b) Multiband applicability of the generalized Szegő limit theorem, and the limitation attributed to the Amoeba formulation.** Verbatim (Abstract of v2; the v1 abstract ends one clause earlier):

> "The Amoeba formulation addresses this challenge by computing the spectral potential rather than the spectrum itself. Based on the (strong) Szegö limit theorem and its topological generalization, this approach reduces the evaluation of the potential to an optimization problem involving the Ronkin function. However, while the generalized Szegö limit theorem is formally applicable in arbitrary dimensions, its implementation is limited to single-band systems, and its applicability to multiband systems remains unclear even in one-dimensional systems."

Verbatim (Section I, v2 HTML; identical wording in v1):

> "While the generalized Szegö limit theorem is, in principle, applicable in arbitrary spatial dimensions, the current Amoeba formalism based on this theorem is limited to single-band systems. In multiband settings, the Ronkin function mixes contributions across bands; when degenerate states with different localization lengths occur, these contributions compete during the optimization of the Ronkin function, causing the standard Amoeba formulation to fail even in one dimension. This breakdown is unavoidable in the presence of symmetry-protected degeneracies, such as for a transpose-type time-reversal symmetry (TRS†), where Kramers pairs localize on opposite edges."

> "In multiband settings, one might hope to resolve the failure of the Amoeba formulation by decomposing the Ronkin function into bandwise contributions. However, a straightforward mathematical decomposition is generally not feasible: the mathematical properties of the Ronkin function are typically not preserved under bandwise separation."

The explicit statement that the criterion is expressed through **Wiener–Hopf partial indices** (verbatim, Section I, v2 HTML):

> "This intimate connection between WHF and topology is central to the generalization of the original theorem to topologically nontrivial Hermitian Hamiltonians, a result known as the Modified Szegö limit theorem [3]."

> "Inspired by these developments, we apply the WHF in this paper to non-Hermitian systems combined with Hermitian doubling [11],[30]. We demonstrate that the WHF of a non-Hermitian non-Bloch Hamiltonian faithfully captures its topological properties. Moreover, this framework allows us to identify the WHF partial indices as the exact criteria determining when the generalized Szegö limit theorem remains valid in multiband systems. By reinterpreting the generalized Szegö limit theorem as a corollary of the modified Szegö theorem, we derive the necessary correction terms to the optimization formalism."

Also, on the WHF/partial-index topological interpretation (verbatim, Section I, v2 HTML):

> "Recently, Ref. [1] demonstrated that the WHF itself offers a direct topological interpretation for one-dimensional Hermitian Hamiltonians: the WHF partial indices count the number of edge modes and thereby establish the bulk-boundary correspondence."

Verbatim (Section II.3, from the ar5iv rendering), the concrete failure mode in class AII†:

> "In the ℤ2-nontrivial phase, two Kramers-related eigenstates localize on opposite edges with inverse localization lengths ±μ1. Consequently, their contributions compete in the minimization in Eq. (14) and the conventional optimization formalism cannot be used."

**(c) Fragments from the later sections (Sections IV–VII) that could only be recovered through the search index** (the prose around them could not be fetched; each snippet is one indexed sentence, as returned by `web_search`):
- "The above discussion illustrates how nonzero winding numbers modify the asymptotic behavior of" — `https://arxiv-org.ezproxy.obspm.fr/html/2511.11349v2#4` (indexed anchor #4; Section III/IV region).
- "It is instructive to interpret the correction term in Eq" — `https://journals.aps.org/prresearch/pdf/10.1103/s43l-h6z6#9#4` (published PRResearch version, page 9 paragraph 4).
- "To this end, we consider a model in which there is a parameter regime, where the partial indices K[σ_μ] d…" — `https://journals.aps.org/prresearch/pdf/10.1103/s43l-h6z6#9#5` (published version, page 9 paragraph 5).
- "Having confirmed that κ serves as the correct topological index, we now verify the central mechanism of our WHF formula…" — `https://browse-export.arxiv.org/pdf/2511.11349#8#6` (arXiv PDF, page 8 paragraph 6).
- "where μ_k and μ_{k+1} denote, respectively, the dominant roots of τ_+^{symp} located just insi…" — `https://browse-export.arxiv.org/pdf/2511.11349#8#5`.
- "where β₁,…,β_{2p} lie inside the unit circle (|β_j| < 1) and β_{2p+1},…,β_{2(p+q)…" — `https://browse-export.arxiv.org/pdf/2511.11349#8#4`.
- "(-1)^{ν[A]} = Pf D(-1)/Pf D(1) = (-1)^{Σ_{κ∈K+[A]} κ}" — `https://browse-export.arxiv.org/pdf/2511.11349#8#3`.
- "To verify the accuracy of the localization length, Fig" — `https://journals.aps.org/prresearch/pdf/10.1103/s43l-h6z6#9#7`.

**Which WSW reference number:** In arXiv:2511.11349 the amoeba references form the group [82]–[87], and "its topological generalization" is cited as **[83]**. The author–year rendering on ar5iv resolves that group as `[Banerjee et al. 2023]` (bib.bib82), `[Wang et al. 2024b]` (bib.bib83), `[Wang 2024]` (bib84), `[Xiong and Hu 2024]` (bib85), `[Xiong et al. 2024]` (bib86), `[Hu 2025]` (bib87). So **WSW = Ref. [83]** (High confidence, mapping via ar5iv author–year labels).

---

### A.2 arXiv:2502.17931 — "Symplectic-Amoeba formulation of the non-Bloch band theory for one-dimensional two-band systems" (Kaneshiro & Peters; published in Phys. Rev. Research)

Sources used: `https://arxiv.org/html/2502.17931` (renders v2, 27 Jun 2025) and `https://ar5iv.labs.arxiv.org/html/2502.17931` (renders an earlier version: it says "band-resolved Ronkin function" where v2 says "symmetry-decomposed"). **Version must be noted, because the reference numbering differs between them.**

**(a) The criterion, including the words "central hole closing" / "hole closing".** Verbatim, Section II.3 (ar5iv rendering):

> "In the Amoeba formulation, the absence of an intermediate region where the Ronkin function is constant indicates that E lies within the spectrum, and the μ*, which minimizes the function, yields the corresponding eigenstate's inverse localization length. This criterion, known as "hole closing," is consistent with the GBZ condition in Eq. (6). The Ronkin function can similarly be defined in higher-dimensional systems. Hence, the Amoeba formulation provides a way to generalize non-Bloch band theory to higher-dimensional systems."

(The arXiv HTML v2 wording of the same criterion is "This criterion, known as "central hole closing," is consistent with the GBZ condition in Eq. (6).")

**(b) How it states the WSW optimization problem and Ronkin function** (verbatim, Section II.3, ar5iv; Eqs. (12)–(13)):

> "Szegö's limit theorem is invalid for topologically nontrivial E, though a generalization is proposed in [Wang et al. 2024b]. According to this conjecture, the potential φ(E) in one-dimensional single-band systems (M = 1) is given by a test potential Φ(E), which is determined through the following optimization problem:
>   Φ(E) = min_μ R_E(μ),   (12)
> where R_E is the Ronkin function defined as
>   R_E(μ) = ∮_{ln|β|=μ} dβ/(2πiβ) ln|ChP(E, β)|,   (13)"

> "This representation indicates that the Ronkin function is determined by the pole and the roots of the characteristic equation enclosed by the circle |β| = e^μ [Xiong and Hu 2023]:
>   R_E(μ) = ln|C_E| − p μ + Σ_{j=1}^{2p} μ_j + Σ_{j=1}^{2p} (μ − μ_j)·θ(μ − μ_j),   (15)
> where μ_j = ln|β_j(E)| and θ(x) is the step function … Thus, the Ronkin function is a convex and piece-wise linear function, and its derivative is quantized to integer values in one-dimensional systems."

> "Using this representation, we can expand the Ronkin function around its minimum as
>   R_E(μ) = a + { −(μ − μ_p) for μ_{p−1} ≤ μ ≤ μ_p ;  0 for μ_p ≤ μ ≤ μ_{p+1} ;  μ − μ_{p+1} for μ_{p+1} ≤ μ ≤ μ_{p+2} },   (16)"

**(c) The limitation attributed to the (WSW) Amoeba formulation — verbatim:**

Abstract (arXiv HTML v2):

> "While this theory provides novel insights into non-Hermitian physics, challenges arise from the multiband nature and symmetry-protected degeneracies, even in one-dimensional cases. In this work, we investigate one-dimensional two-band class AII† systems, where Kramers pairs invalidate the conventional Amoeba formalism."

Section I (ar5iv):

> "The Amoeba formulation defines the "electrical" potential of the spectrum, making it possible to access the density of states (DOS) under OBC for single-band class A systems in arbitrary dimensions. In particular, the generalized Szegö's limit theorem reduces the potential computation to an optimization of the Ronkin function, which is obtained by the PBC spectrum. Furthermore, as the Ronkin function takes the inverse localization length as an argument, the solution to this optimization problem returns the correct localization length."

> "While the Amoeba formulation provides a powerful framework for developing the non-Bloch band theory of one-band systems, generalizing it to multiband systems remains challenging. In particular, degeneracies arising from the transpose-type time-reversal symmetry (TRS†) impose constraints on the Ronkin function [Wang 2024], which obstruct both the optimization problem and the application of the generalized Szegö's limit theorem. These constraints are especially significant in symplectic classes, where TRS† leads to Kramer's degeneracies, further complicating the Amoeba formulation. Consequently, even for one-dimensional systems, the Amoeba formulation fails to apply to symplectic classes."

Section III, opening sentence (this is the sentence the coordinating agent saw as "…discuss the limitations of the conve…"):

> "In this section, we focus on one-dimensional two-band class AII† systems and discuss the limitations of the conventional Amoeba formalism."

Section III.2, "Breakdown of the Amoeba formulation …" (ar5iv; this is the sharpest stated failure of the WSW optimization in a multiband case):

> "Since the Ronkin function becomes even [Wang 2024], μ* = 0 always minimizes the Ronkin function. This seemingly implies that all bulk modes are delocalized, and the NHSE is absent in class-AII† systems. However, this conclusion contradicts both theoretical analyses and numerical results, which consistently demonstrate the existence of localized bulk modes and the NHSE in such systems …"

> "For class AII† systems, we can calculate the Ronkin function in Eq. (13) as
>   R_E(μ) = ln|C_E| − 2pμ + Σ_{j=1}^{2p}(μ + μ_j)·θ(μ + μ_j) + Σ_{j=1}^{2p}(μ − μ_j)·θ(μ − μ_j),   (24)
> where μ_j = ln|β_j(E)|. Note that all μ_j are positive as required by Eq. (22)."

**(d) Which WSW reference number.** In the arXiv HTML (v2) numbering, the amoeba group is [58]–[63] and the phrase "a generalization is proposed in [59]" plus "the Amoeba formulation [59] in class A systems" identifies **WSW = Ref. [59]**. In the ar5iv rendering (older version numbering), the same citations appear as `[Wang et al. 2024b](#bib.bib58)` — i.e. **WSW = Ref. [58]** there. Report whichever version's numbering you are using. (Also cited in that paper: "the Ronkin function is determined by the pole and the roots … [61]" in v2 numbering, `[Xiong and Hu 2023]` = bib.bib60 in the ar5iv numbering — again a one-off version shift.)

---

## B. Xiong, Xing & Hu — arXiv:2407.01296, "Non-Hermitian skin effect in arbitrary dimensions: non-Bloch band theory and classification" (published as "Geometry-adaptive formulation of non-Bloch bands in arbitrary dimensions and spectral instability", Commun. Phys. **9** (2026) 2546, DOI 10.1038/s42005-026-02546-2)

Sources used: `https://arxiv.org/html/2407.01296`, `https://ar5iv.labs.arxiv.org/html/2407.01296`. The Nature/Springer versions could not be fetched (see route notes above). Section V itself lies beyond the truncation point; its content below is quoted from the search index where marked.

**(a) The critique of WSW — verbatim (Section I, ar5iv rendering, in which the reference appears as "Wang et al. 2024" = bib.bib65):**

> "A recent work proposed a higher-dimensional generalization of the GBZ condition via the Amoeba formulation [Wang et al. 2024], which neglects geometric information and yields geometry-irrelevant non-Bloch spectra. Yet, insights from 1D NHSE suggest that boundary conditions or lattice geometry play a crucial role in the asymptotic spectral structure and skin modes. As will be demonstrated in this paper, in higher-dimensional non-Hermitian systems with clear boundaries, incorporating geometric information is essential to formulating the non-Bloch band theory."

> "Using Ronkin's function, the GBZ is proposed [Wang et al. 2024] to be a dD object embedded in a 2dD space for dD non-Hermitian systems without specifying the geometric information."

> "In Ref. [Zhang et al. 2022b], varying spectral densities associated with different lattice shapes are reported. However, numerical results in Ref. [Wang et al. 2024] suggest that systems under weak perturbations or with smooth boundaries should have the same DOS in the TDL. Refs. [Hu 2023]; [Wang et al. 2024]; [Xu et al. 2024] further conjecture that the energy spectra for any lattice geometry in the TDL should be given by the spectra from the Amoeba formulation. This debate seems difficult to settle since numerical errors arising from the non-normality of the Hamiltonians are inevitable."

(Note: the "neglects geometric information" clause is the critique the coordinating agent had seen; it is in the **Introduction**, not in Section V.)

**(b) Their preview of Section V — verbatim (Section I, ar5iv):**

> "Section V provides a pedagogical introduction to the Amoeba formulation. Through an analytically tractable model, we further benchmark our formalism and show its disparities with Amoeba. The relationships between different types of energy spectra and their geometric dependencies are established through a spectral theorem and the first conjecture."

**(c) The d-dimensional Ronkin/amoeba statement, and how they restate the Szegö/Ronkin chain — the 1D version is inside the fetched range (verbatim, Section II.2, ar5iv):**

> "To obtain the non-Bloch spectra, the key observation is that in the TDL, the potential function φ(E) takes the integral form:
>   lim_{N→∞} φ(E) = ∫₀^{2π} dk/2π log|det[H(e^{ik+μ}) − E]|,   (12)
> according to Szegö limit theorem. Here, the parameter μ should be properly chosen such that |β_p(E)| < e^μ < |β_{p+1}(E)|. Notably, this condition is not satisfied if E is located on the spectral arcs due to |β_p(E)| = |β_{p+1}(E)|."

> "… one can obtain the local form of the electrostatic potential [Xiong and Hu 2024]:
>   lim_{N→∞} φ(E) = Σ_{j=p+1}^{p+q} log|β_j(E)| + log|t_q|,   (13)
> where β_j are the zeros of ChP sorted as |β₁| ≤ |β₂| ≤ ⋯ ≤ |β_{p+q}|. Thus the potential φ(E) has contributions from the q roots of the largest moduli."

> "By scanning all possible deformation parameters μ ∈ (−∞, +∞), it can be rigorously proven … that
>   lim_{N→∞} φ(E) = min_μ ∫₀^{2π} dk/2π log|det[H(e^{ik+μ}) − E]|, ∀E.   (14)
> This means that the local potential at E is the minimum among all possible deformed spectral loops. Thus we propose the following principle of non-Bloch bands from the perspective of potential landscape.
> Guiding principle: Among all possible spectral deformations, the one with the minimum spectral potential corresponds to the potential generated by the non-Bloch bands."

> "The GBZ condition is derived from these boundary constraints. However, in d ≥ 2D, the number of boundary equations is of O(L^{d−1}) order (L is the linear length of the lattice) [See Fig. 2(a2)], which tends to infinity in the TDL. Thus, it is impossible to simultaneously solve these boundary equations and determine the GBZ condition."

**(d) Section V's d-dimensional amoeba formula — recovered only as an indexed snippet (verbatim, single sentence/equation as returned by `web_search`):**

> "Φ_Amoeba(E) = min_𝝁 ∫ d^d𝒌/(2π)^d log|f(e^{i𝒌+𝝁}, E)|"

Source URL as reported by the index: `https://arxiv-org.ezproxy.obspm.fr/html/2407.01296v1#5` (the index's anchor #5; the corresponding section is the one the paper calls its pedagogical introduction to the Amoeba formulation). **I could not fetch the surrounding prose of Section V**, so this is the only WSW-formalism equation of theirs I can quote verbatim, together with:

> "In our formulation, the geometric information is input through basis transformation, and the spectral potential under a specific…"

from `https://ar5iv.labs.arxiv.org/html/2407.01296#5`.

The published version's Supplemental Information is indexed with the companion 2D expression (fragment only):

> "Φ_Amoeba(E) = ∬ dk₁dk₂/(2π)² log |f(e^{ik₁+μ₁,min}, e^{ik₂+μ₂,…"

from `https://static-content.springer.com/esm/art%3A10.1038%2Fs42005-026-02546-2/MediaObjects/42005_2026_2546_MOESM2_ESM.pdf#3#2` (a PDF I cannot fetch; the fragment comes from the search index only).

**(e) Which WSW reference number:** **Ref. [65]** — confirmed independently of the numbered text: in the ar5iv rendering the same sentence's citation is the hyperlink `[Wang et al. 2024](#bib.bib65)`, and the numbered arXiv HTML/PDF text of the same sentence reads "[65]". (High confidence; this matches the coordinating agent's note.)

---

## C. Yang & Fang — arXiv:2503.11505, "Real-time edge dynamics of non-Hermitian lattices" (Phys. Rev. Lett. **135**, 186401 (2025))

Source used: `https://arxiv.org/html/2503.11505v2` (fetch truncated inside the End Matter / Appendix E).

**Finding: no restatement or extension of the amoeba/Ronkin formalism appears in the readable portion.** The words "amoeba" and "Ronkin" do not occur anywhere in the text I was able to fetch (abstract, Introduction, "The saddle point method", "Edge dynamics", and the first part of the appendices/supplemental material). The paper's higher-dimensional discussion is instead framed as saddle-point/Lefschetz-thimble theory. Its only statement about the state of higher-dimensional non-Bloch theory is (verbatim, Introduction):

> "While analytic methods exist, they are computationally demanding [47], and generalization to higher-dimensional cases remains an ongoing area of research [27], [51], [52], [53], [54]."

I could **not** read this paper's reference list (it lies beyond the fetch truncation point), so I cannot confirm or deny whether WSW is among refs [51]–[54]. Verbatim, its own contrast of the eigenstate-based vs. dynamical description (relevant to any "what the amoeba formulation does not capture" argument):

> "Notably, the DSP z_s need not be on the GBZ or in the point gap, and the stationary state is not a skin mode, defying previous expectations [56], [57], [58], [59]."

> "These results demonstrate that for non-Hermitian systems, real-time dynamics can display qualitatively different physics from what eigenstate analysis predicts [68]."

and its explicit caveat on the topological equivalence of BZ and GBZ (a stated unproven assumption):

> "More precisely, "topologically equivalent" means that the BZ and GBZ are homotopic, i.e., can be continuously deformed into each other. This fact is widely believed to hold. We are not aware of any proofs, nor of any counterexamples."

---

## D. Two further derivative papers

### D.1 arXiv:2608.28577 — Gu, Fu, Hu & Wang, "How Long-Range Tails Reshape Non-Hermitian Spectra" (Zhong Wang is a WSW author)

Source used: `https://arxiv.org/html/2608.28577v1` (a short Letter + Supplemental Material; the amoeba review section is fully inside the fetched range). Identity confirmed at `https://arxiv.org/abs/2608.28577`.

Verbatim (abstract):

> "In two or higher dimensions, we formulate a squeezed amoeba formulation describing the reconstructed spectral density of states (DOS)."

Verbatim ("DOS in 2D from squeezed-domain Ronkin function"):

> "It is naturally expected that the competing mechanism also emerges in higher-dimensional systems with small long-range hoppings. However, the higher-dimensional GBZ cannot be calculated directly. Instead, the OBC spectrum is encoded in an amoeba formulation of higher-dimensional non-Bloch band theory [3]. Therefore, a generalized amoeba theory compatible with long-range perturbations is needed."

> "To this end, we first review the standard amoeba formulation for a short-range 2D non-Hermitian Hamiltonian H₀ described by a non-Bloch Hamiltonian h₀(β_x, β_y)."

> "A central concept in the amoeba formulation is the Ronkin function, defined as a 2D integral over a torus parametrized by real variables θ_{x,y} ∈ [0, 2π]:
>   R_E(μ_x, μ_y) = ∫_{T²} dθ_x dθ_y/(2π)² log|f(e^{μ_x+iθ_x}, e^{μ_y+iθ_y})|,   (5)
> where f(β_x, β_y) = det[E − h₀(β_x, β_y)] is the charateristic equation for a given energy E."

> "In the amoeba formulation, φ(E) corresponds to the global minimum of the Ronkin function:
>   φ(E) = min_{μ_x, μ_y ∈ ℝ} R_E(μ_x, μ_y).   (6)"

> "Ref. [3] has shown that R_E(μ_x, μ_y) is a convex function over the entire (μ_x, μ_y) plane, whose minimum is attained either at a single point or on a finite plateau. It has been established that if the Ronkin function minimum forms a plateau, the corresponding energy E does not belong to the OBC spectrum of H₀. However, if the minimum is attained at a single point (μ̃_x, μ̃_y), the spectral DOS at E is nonzero and the corresponding eigenstate exhibits an asymptotic form |ψ_E(x, y)| ∼ e^{μ̃_x x + μ̃_y y}. This correspondence implies that any modification of the Ronkin-function minima directly affects the spectral DOS ρ(E)."

**Which WSW reference number:** unverified. In this paper the amoeba formulation is cited as "[3]" but the link anchors point into the **Supplemental Material bibliography** (`#as1_bib.bib3`), and the supplemental reference list is beyond the fetch truncation point. So "the amoeba formulation [3]" cannot be safely equated with WSW from what I could read; it may be Banerjee et al. or WSW. Flag for the coordinator.

### D.2 arXiv:2607.22976 — Xu, Wang, Deng & Yi, "Spectral Topology and Non-Bloch Band Theory for Domain-Wall Systems"

Source used: `https://arxiv.org/html/2607.22976v1`. Identity confirmed at `https://arxiv.org/abs/2607.22976`. This paper **extends** the Ronkin-function formalism to a domain-wall ring; it is a direct descendant of the WSW machinery and states the WSW-type single-domain objects as follows.

Verbatim (abstract):

> "We then obtain the conditions for the generalized Brillouin zone (GBZ) in the complex momentum space, by extending the Ronkin-function formalism to the domain-wall configuration."

Verbatim (main text):

> "Further analysis through the Ronkin functions in the complex momentum space [6], [17], [18], [33], [34], [36], [35] reveals that the DW spectrum is generally divided into the standing-wave- and traveling-wave-like sectors …"

> "Inspired by the single-domain Ronkin-function formalism under the OBC [17], we characterize this critical limit by the constrained Ronkin function."

> "We define the constrained Ronkin function for an n-DW ring
>   ℛ(𝝁; E) = Σ_{α=1}^{n−1} r_α R_α(μ_α; E) + r_n R_n(−Σ_{α=1}^{n−1} (r_α/r_n) μ_α; E),   (8)
> where the Ronkin function in domain α is
>   R_α(μ_α, E) = ∫₀^{2π} dk_α/2π log|f_α(β_α, E)|,   (9)
> and β_α := e^{ik_α+μ_α}. For each domain, f_α(β_α, E) = det[h_α(β_α) − E] gives a Laurent series in β_α …"

> "By construction, R_α(μ_α, E) is a piecewise linear function, and its derivative with respect to μ_α equals the spectral winding of the corresponding non-Bloch Hamiltonian [17]. Hence, we have ∂_{μ_α} R_α(μ_α, E) = w_α(μ_α, E) …"

> "Importantly, the derivative of the constrained Ronkin function Eq. (8) with respect to μ_α is (α ≠ n)
>   ∂_{μ_α} ℛ = r_α [ w_α(μ_α, E) − w_n(μ_n, E) ],   (11)
> so that the solutions to Eq. (7) corresponds to a flat region in ℛ(𝝁, E), characterized by ∂_{μ_α} ℛ = 0, within an (n−1)-dimensional 𝝁 space (that is, the solution space)."

> "Typically, a two-dimensional flat region of the Ronkin function emerges, as in (b), indicating that E₁ lies in a point gap of the DW-ring. The collapse of the flat region, illustrated in (c) and (d), places the corresponding E₂ and E₃ on the DW eigenspectrum."

**Which WSW reference number:** unverified. The single-domain Ronkin formalism is cited as "[17]" (plain bracketed number); the bibliography lies beyond the truncation point, and no author–year rendering was available for this paper (it has no ar5iv render). "[17]" is plausibly WSW or Xiong–Hu; I cannot distinguish them from what I read.

---

## E. Appendix — WSW's own later-section fragments, and one further census fragment

These are **not** derivative papers; they are included only because they were surfaced by the search index during this session and corroborate WSW's own wording/equation environment. Each is a verbatim indexed fragment with its source URL; I did not read the surrounding prose.

From the WSW PRX PDF (`https://journals.aps.org/prx/pdf/10.1103/PhysRevX.14.021011`, page 9 anchors):
- "With the results of the DOS in the previous section, the proof of this theorem is now simple" (`#9#5`)
- "We now apply the explicit formula of the Ronkin function to \operatorname*{det}[E − h(β)] = a_{−M}(E)β^{−M} + \dots…" (`#9#3`)
- "Taking a circle |E| = R that surrounds the whole PBC spectrum, the average potential on this circle is \tilde{\Phi}_{\math…" (`#9#7`)
- "This model is known to have a Chern insulator phase as well as a trivial insulator phase in the Hermitian limit (γ = 0)…" (`#9#6`)

From WSW's arXiv HTML rendering (`https://ar5iv.labs.arxiv.org/html/2212.11743v3#3`):
- "where is a holomorphic function with gg, and g(0) ≠ 0 (z_k) are the zeros of k = 1,⋯,l enclosed by the cir…" (index snippet, LaTeX-degraded)

From a further derivative paper that uses the amoeba-hole vocabulary (`Phys. Rev. Research` **7**, 023233, "Anatomy of higher-order non-Hermitian skin and boundary modes", `https://journals.aps.org/prresearch/pdf/10.1103/PhysRevResearch.7.023233#14#8`):
- "where \(\mathcal{V}_{\mathrm{hole}}\) denotes the volume of amoeba hole" (fragment only; the paper's use of WSW's hole variable could not be read further).

---

## F. Citation-number corroboration summary

| Derivative paper | WSW cited as | Evidence (what I actually saw) | Confidence |
|---|---|---|---|
| Xiong, Xing & Hu, arXiv:2407.01296v1 | **Ref. [65]** | Numbered text "[65]" for the Amoeba formulation; ar5iv inline link for the same sentence resolves to `[Wang et al. 2024](#bib.bib65)` | High |
| Kaneshiro & Peters, arXiv:2511.11349 (v1 & v2) | **Ref. [83]** | Amoeba group [82]–[87]; "its topological generalization [83]"; ar5iv resolves bib.bib83 = "Wang et al. 2024b" | High |
| Kaneshiro & Peters, arXiv:2502.17931 | **[59]** in `arxiv.org/html` (v2) numbering; **[58]** in the ar5iv rendering | v2: "the Amoeba formulation [59]", "a generalization is proposed in [59]"; ar5iv: `[Wang et al. 2024b](#bib.bib58)` | Medium–High (version-dependent — always state the version) |
| Gu, Fu, Hu & Wang, arXiv:2608.28577 | "[3]" pointing into the **supplemental** bibliography (`#as1_bib.bib3`) | Anchor target only; bibliography not readable | Unverified |
| Xu, Wang, Deng & Yi, arXiv:2607.22976 | "[17]" (single-domain Ronkin function under OBC) | Plain bracketed number; bibliography not readable | Unverified |
| Yang & Fang, arXiv:2503.11505 | not found | No occurrence of "amoeba"/"Ronkin" in the readable portion | — |

---

## G. What I could NOT read (explicit list)

1. **Sections III–VII of arXiv:2511.11349 (Kaneshiro & Peters)** — i.e. the WHF section, the Strong-Szegő section, and the central Section V where the multiband applicability criterion in terms of partial indices is stated, together with the "correction terms". Their HTML text is ~3× the fetch budget because of MathML duplication; PDF routes are unsupported. Only indexed one-sentence fragments (listed in A.1(c)) are available, and none of them contains the complete criterion sentence. **This is the single largest gap in this file.**
2. **Section V ("Amoeba formulation") of arXiv:2407.01296 (Xiong, Xing & Hu)**, and its subsections "Amoeba spectra vs non-Bloch spectra" and "Spectral relations": the prose, the hole-closing discussion, and the equations other than the single `Φ_Amoeba(E) = min_𝝁 ∫ d^d𝒌/(2π)^d log|f(e^{i𝒌+𝝁},E)|` snippet. Also unread: their **two conjectures** as literally worded, and the published Communications Physics version (nature.com redirected to `idp.nature.com`; link.springer.com unreachable; not in Europe PMC by DOI).
3. **arXiv:2502.17931 (Kaneshiro & Peters) Sections III.2–V** beyond the first paragraph of III.2 that the ar5iv fetch reached — i.e. the band-resolved/symmetry-decomposed Ronkin construction in full, the Legendre-transformation step, and their restated generalized Szegő limit theorem as an equation. Also unread: Appendix A (why band-resolved Ronkin functions fail to inherit the mathematical properties of the total one) and Appendix B (GBZ criterion from the spectral potential).
4. **arXiv:2608.28577 (Gu et al.) Supplemental Material** beyond its opening: in particular the section "Ronkin function formulation of squeezed GBZ in 1D" was cut off mid-equation, and the **supplemental reference list** (needed to identify "[3]") was never reached.
5. **arXiv:2607.22976 (Xu et al.) Sections III.3–VI** and its reference list: the travelling-wave GBZ conditions (Eqs. (12)–(14) region, cut off at "Case I, ∃w, such that…"), the spectral-winding section, and the bibliography (needed to identify "[17]").
6. **Reference lists in general**: for every paper above, the bibliography sits at the end of the HTML and is therefore beyond the truncation point. Every WSW reference number reported in section F was reconstructed from (i) numbered in-text citations and (ii) ar5iv author–year hyperlink targets — **never from a bibliography I actually read.** The two "Unverified" rows are cases where neither route was available.
7. **Abstracts/full texts of the published journal versions**: APS PDFs (`unsupported content type`), APS harvest fulltext endpoint (also PDF; `5s1z-5t9r` additionally failed at the transport level), Nature HTML (auth redirect), and Springer ESM PDFs (unsupported). Page anchors such as `#9#4` refer to the PDFs I could not open.
8. **Yang & Fang's "Notes added" section and full appendices** (fetch truncated inside Appendix E); consequently I cannot state with certainty that the published PRL never mentions the amoeba formulation — only that it does not appear in the ~50 000 characters I could read, which covered the abstract, Introduction, "Saddle point method", "Edge dynamics", and the beginning of the supplemental material.

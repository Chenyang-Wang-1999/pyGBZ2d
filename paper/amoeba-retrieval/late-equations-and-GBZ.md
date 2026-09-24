# WSW (arXiv:2212.11743v3 / PRX 14, 021011) — Late equations, appendix eqs. (a11)/(a12), Szegő/Widom machinery, asymptotics and GBZ

**Target:** Hong-Yi Wang, Fei Song, Zhong Wang, "Amoeba Formulation of Non-Bloch Band Theory in Arbitrary Dimensions", arXiv:2212.11743v3, Phys. Rev. X **14**, 021011 (2024).

**Retrieval date/context:** delegated retrieval task. Methods used: `web_search` phrase mining (the search index holds the publisher PDF and the arXiv/ar5iv HTML with page/paragraph anchors), `web_fetch` of the APS abstract page (its reference/footnote block renders *verbatim* article text), `web_fetch` of arXiv HTML / ar5iv HTML / APS abstract page, and `web_fetch` of derivative papers.

**Honesty preamble (read this first).** `web_fetch` returns only the first ~50,000 characters of any document. Direct fetches of the article body (arXiv HTML `arxiv.org/html/2212.11743v3`, ar5iv `ar5iv.labs.arxiv.org/html/2212.11743v3`) were **truncated inside Section III / Section II.2** on every attempt, so Sections IV.3, V, VI, VIII and Appendix A are **not readable by direct fetch**. All WSW body text below therefore comes from either (i) the section of the paper that *is* within the first 50 KB, or (ii) verbatim passages surfaced as titles of search-index entries, or (iii) the footnote/endnote block of the APS abstract page. I mark the provenance of every single quote. Nothing is paraphrased as if it were a quote; where I only have a fragment, I say so.

---

## PART I — WSW's OWN WORDS (verified verbatim)

### I.1 Footnote / endnote block of the APS abstract page — FULLY READABLE

Source: https://journals.aps.org/prx/abstract/10.1103/PhysRevX.14.021011 (fetched; the numbered reference list on that page includes the article's footnotes, rendered verbatim). The page supplies anchors `#d44` (Eq. (44)), `#da11`, `#da12`.

I.1.a — **Endnote 47 (recovered verbatim, with its equation links).** Note the text states the symbol decomposition inline, not as a display equation:

> "In the presence of EPs, we have an additional nilpotent term N in the spectral decomposition: H=∑_n ε_n|n_R⟩⟨n_L|+N. The Green's function also acquires an additional term, roughly as (E−H)^{-1}=∑_n|n_R⟩⟨n_L|/(E−ε_n)+∑(nilpotent matrix)/(E−ε_n)^p, where p≥2 [1]. The final term vanishes after the contour integral in Eq. (44), and therefore, Eq. (44) still holds."

(Wording differences from the task brief: the published text reads "In the presence of EPs" and "The final term vanishes after the contour integral in Eq. (44)"; it does **not** contain the phrase "vanishes after the contour integral" as a standalone string, and does not say "we have an additional nilpotent term N" without "In the presence of EPs". Reference "[1]" is Ashida–Gong–Ueda, Adv. Phys. **69**, 249 (2020).) ⇒ **Eq. (44) is a contour integral of the Green's function that survives the presence of exceptional points.**

I.1.b — **Endnote 61 (recovered verbatim in full, with its equation links `#da11`, `#da12`):**

> "As a heuristic example, consider σ(e^{iθ})=ε+e^{iθ}. For ε=10, the symbol is homotopic to a constant symbol, and the winding number is zero. One can check that both Eqs. (a11) and (a12) hold. As a comparison, the ε=1/10 case has winding number 1, meaning that the symbol is not homotopic to a constant symbol. Accordingly, Eq. (a11) is satisfied, while Eq. (a12) is not. Instead, it is homotopic to a unitary translation operator, whose symbol is σ(e^{iθ})=e^{iθ}. It is now obvious that a unitary translation maps the leftmost mode to zero, thus having one zero singular value. Being in the same topological sector, the matrix with ε=1/10 also has a zero singular value (with a correction exponentially small in size)."

⇒ **Eqs. (a11) and (a12) are in an Appendix** (typeface "a11"/"a12" = appendix equations; the only appendix of the paper is "APPENDIX A: A brief proof of Szegő's limit theorem"). The brief's phrasing "maps the leftmost mode to zero and thus has one zero singular value" is close but not identical to the published text ("maps the leftmost mode to zero, thus having one zero singular value").

I.1.c — **Endnote 34 (recovered verbatim):**

> "We also check that adding local randomness in the bulk has a similar effect. Here, the bulk randomness should be sufficiently weak such that the resulting energy spectrum reflects the properties of the pristine Hamiltonian. Specifically, one can add random potential only to N_d unit cells, with N_d/N→0 (though N_d→∞) when taking the large-size limit N→∞ (N is the number of unit cells)."

I.1.d — **Endnote 38 (recovered verbatim):** "When taking the zero-coupling limit, one sets the chain length to be large and fixed."

I.1.e — **Endnote 46 (recovered verbatim):** "Note that possible boundary terms are omitted, in the same spirit as in 1D."

I.1.f — **Attribution caveat on the sentence in I.1.e.** In the page's rendered list, the entry numbered **46** carries the string "Note that possible boundary terms are omitted, in the same spirit as in 1D." I did not see any *other* endnote carrying that same sentence, and I could not independently verify a separate endnote **44** with this wording. **Net position:** the sentence is verbatim WSW, and the rendered list assigns it to position 46 (consistent with the task brief); the possibility that it is 44 rather than 46 could not be fully excluded from the rendered payload, so treat the number as **46 per the rendered list**, not as independently confirmed.

### I.2 Article body text readable within the fetch window

Source: https://arxiv.org/html/2212.11743v3 (fetched; truncated inside Sec. II.2) and https://ar5iv.labs.arxiv.org/html/2212.11743v3 (fetched; truncated at the end of Sec. III). These two give **identical** text for the overlap, so the quotes below are doubly confirmed across two independent renderings.

- Section structure confirmed from the arXiv HTML table of contents: "IV Energy Spectra and density of states" → **"IV.1 Statement of the proposal"; "IV.2 Numerical evidence"; "IV.3 Derivation"**; "V Amoeba hole closing and spectral boundary"; "VI Generalized Brillouin zone"; "VII Non-Bloch band topology"; "VIII Spectral inequalities"; "IX Concluding remarks"; "**A A brief proof of Szegő's limit theorem**"; "B Invariance under change of basis".
- The 1D GBZ equation of the paper is **Eq. (9)**: |β_M(E)|=|β_{M+1}(E)|.
- 2D characteristic equation is **Eq. (13)**: det[E−h(β)]=0.
- The 2D single-band model is **Eq. (12)**; the non-Hermitian SSH model is **Eq. (11)**.
- WSW's own words on the hole↔OBC-spectrum heuristic (from Sec. II.2, verbatim): "From the above examples, we observe that the absence (presence) of a hole in the amoeba of the characteristic polynomial could be an indicator of the energy E being (not being) in the OBC energy spectrum. This is a key observation of the present work."
- WSW's own words introducing the amoeba (Sec. II.2, verbatim): "…we plot the solutions to the 2D characteristic equation, followed by mapping (μ_x,μ_y)=(log|β_x|,log|β_y|) … This geometrical object is known as the amoeba in mathematics literature…"

### I.3 Verbatim WSW passages recovered as search-index entry titles (publisher PDF and arXiv PDF)

These are exact in-article strings returned by `web_search` as the *titles* of source entries. Their parenthetical anchors are the index's own page/paragraph markers. **Each is a fragment**: the index shows only the leading portion of the paragraph. I quote only what was returned, and say where it is cut.

Source A: `https://journals.aps.org/prx/pdf/10.1103/PhysRevX.14.021011` (anchors `#9#N`), Source B: `http://export-test.arxiv.org/pdf/2212.11743` (anchors `#9#N`), Source C: `https://arxiv.org/pdf/2503.11505`, Source D: `https://arxiv-org.ezproxy.obspm.fr/html/2212.11743v3` (anchors `#N`).

1. **[Szegő/Widom — the key sentence]** (Source B, `#9#4`; and Source A, `#9#4`): 
   > "Szego's limit theorem was originally established for 1D Hermitian Toeplitz matrices [41], but thereafter generalized by Widom _e..." 
   
   The string is **truncated by the index mid-word ("_e...")**. Reference [41] is G. Szegő, Math. Ann. **76**, 490 (1915); the Widom references in WSW's list are [42] Adv. Math. **13**, 284 (1974) and [43] Adv. Math. **21**, 1 (1976), plus [44] J. Funct. Anal. **39**, 182 (1980) ("Szego's limit theorem: The higher-dimensional matrix case") and [45] Doktorskii, Sib. Math. J. **25**, 701 (1984). The sentence immediately preceding this one is (Source A, `#9#4`, verbatim): 
   > "The language of Toeplitz matrices is very useful in addressing tight-binding Hamiltonians".

2. **[Use of the Ronkin formula on the characteristic polynomial]** (Source A, `#9#3`, verbatim): 
   > "We now apply the explicit formula of the Ronkin function to \(\operatorname *{det}[E - h(\beta)] = a_{- M}(E)\beta^{- M} + \dots..."
   
   Cut off. This confirms WSW feed the *characteristic polynomial's* Laurent expansion a_{−M}(E)β^{−M}+…+a_N(E)β^N into an explicit Ronkin formula — i.e. this is the derivation step behind φ(E)=min_μ R_E(μ).

3. **[Consequence of the Ronkin formula / hole condition]** (Source A, `#9#8`, verbatim): "Therefore, a condition for Eq" — cut off immediately after "Eq".

4. **[Proof of the spectrum theorem from the DOS]** (Source A, `#9#5`, verbatim): 
   > "With the results of the DOS in the previous section, the proof of this theorem is now simple"
   
   This is the opening of the proof passage in Sec. IV.3 (or Sec. V); it confirms the logical order WSW use: DOS (Eq. (44)) first, then the spectral theorem.

5. **[Averaging the potential on a large circle]** (Source A, `#9#7`, verbatim): 
   > "Taking a circle \(|E| = R\) that surrounds the whole PBC spectrum, the average potential on this circle is \(\tilde{\Phi}_{\math..."
   
   Cut off mid-symbol. This is WSW's device for relating the potential far outside the spectrum to the log|leading coefficient| — the same device as the a_{−M}/a_N expansions in item 2.

6. **[Alternative proof of a spectral inequality]** (Source B, `#9#7`, verbatim): "An alternative proof of the spectral inequality is based on Eq" — cut off after "Eq".

7. **[Explicit Szegő-type evaluation in the appendix]** (Source B, `#9#3`): a display formula returned verbatim as 
   `\[R_{g}(\log|R|) =\log|g(0)|+\log\frac{|z_{2}|}{|z_{1}|}+2\log\frac{|z_{3}|}{|z_{2}| }+3\log\frac{|z_{4}|}{|z_{3}|}\] \[+\cdots+...\]`
   
   and the **same identity in a second, independent rendering** (Source D, `#3`, verbatim): 
   > "Amoeba Formulation of Non-Bloch Band Theory in Arbitrary Dimensions - 12π∫02πdθlog|g(Reiθ)|=log|g(0)|+∑k=1llog|Rzk|," 
   
   i.e. R_g(log|R|) = (1/2π)∫_0^{2π} dθ log|g(Re^{iθ})| = log|g(0)| + Σ_{k=1}^{l} log|R/z_k|, with the accompanying statement (Source D, `#3` / ar5iv `#3`, verbatim): 
   > "where is a holomorphic function with gg, and g(0)≠0g(0)\neq 0 (zkz_{k}) are the zeros of k=1,⋯,lk=1,\cdots,l enclosed by the cir..." (index text is garbled by whitespace/math reflow: intended reading is "where g is a holomorphic function with g(0)≠0 and z_k (k=1,⋯,l) are the zeros of g enclosed by the circle").
   
   The *second* form — `log|g(0)| + log|z2/z1| + 2log|z3/z2| + 3log|z4/z3| + ⋯` — shows WSW (or the index's rendering of their appendix) group the zeros **by multiplicity** (coefficients 1, 1, 2, 3, … are the multiplicities), i.e. the Ronkin function of a single-variable polynomial evaluated on a circle equals log|g(0)| plus a multiplicity-weighted telescoping sum of log|R/z_k|. This is exactly the "explicit formula of the Ronkin function" appealed to in item 2. **Caveat:** this display sits in Appendix A (or in Sec. III's discussion of the Ronkin function), and I could not read the surrounding sentences, so I cannot state which equation number it carries, nor whether the exponent-like leading coefficients (1,1,2,3) are multiplicities or are the index's own transcription of β-powers.

8. **[Lattice coordinates in the spectral-inequality section]** (Source B, `#9#6`, verbatim): 
   > "where \(\bm{x}\) are 2D integer coordinates, and \(\bm{e}_{j}\) is the unit vector in the \(j\)th direction [see Fig" — cut off.

9. **[The model whose phase boundary is fixed analytically]** (Source D, `#7`, verbatim): 
   > "For the present model, the precise phase boundary can even be analytically determined by the amoeba formulation" — cut off.

10. **[GBZ derivation uses a point fixed in the amoeba hole]** (Source D, `#5`, verbatim): 
    > "where with 𝜷=e^{𝝁+i𝒌} fixed in the central hole of the amoeba of 𝝁" 
    
    Cut off / reflowed. This is from the GBZ section (or the DOS derivation): WSW fix β = e^{μ+ik} with **μ inside the central hole of the amoeba**, which is the natural way to define the GBZ as a subset of the amoeba's boundary.

11. **[Sec. IV.2 statement about numerical evidence]** (Source D, `#6`, verbatim, fragment): 
    > "Now we provide evidence for Eqs"
    
    Cut off immediately (this is the opening of Sec. IV.2 "Numerical evidence").

### I.4 WSW's own "Popular Summary" (verbatim, APS abstract page)

> "The key mathematical tool in our work is an elegant geometrical object called an amoeba. The amoeba is essentially a shadow or a projection of a polynomial of complex variables onto real space. This name is suggested by its typical appearance, featuring elongated "tentacles" and sometimes "vacuoles" in its body. We show that a surprisingly rich amount of information about a non-Hermitian system is encoded in an associated amoeba. For instance, the energy spectrum can be straightforwardly obtained by looking at whether a certain vacuole appears inside the amoeba."

### I.5 WSW's own Introduction (verbatim, arXiv HTML, read directly)

> "…despite the geometry-dependent NHSE, the existence of a universal spectrum (amoebic spectrum) to which the OBC spectrum under any generic geometry converges. Furthermore, the amoeba inspires a definition and the associated algorithm of the GBZ in arbitrary spatial dimensions."

---

## PART II — DERIVATIVE-PAPER RESTATEMENTS (clearly NOT WSW's own words)

### II.1 Gu–Fu–Hu–Wang, arXiv:2608.28577 ("How Long-Range Tails Reshape Non-Hermitian Spectra"), Zhong Wang is a co-author

Source (fetched in full past the truncation point, because the Letter is short): https://browse-export.arxiv.org/html/2608.28577v1

Verbatim, their Sec. "DOS in 2D from squeezed-domain Ronkin function":

> "A central concept in the amoeba formulation is the Ronkin function, defined as a 2D integral over a torus parametrized by real variables θ_{x,y}∈[0,2π]:
> R_E(μ_x,μ_y)=∫_{T²} (dθ_x dθ_y)/(2π)² log|f(e^{μ_x+iθ_x}, e^{μ_y+iθ_y})|,   (5)
> where f(β_x,β_y)=det[E−h_0(β_x,β_y)] is the charateristic equation for a given energy E."

> "In the thermodynamic limit with generic (i.e., non-fine-tuned) open boundary conditions, the spectral density of states (DOS) ρ(E) is given by ρ(E)=1/(2π) ∇²φ(E), where a potential function φ(E) can be viewed as a 2D Coulomb potential on the complex energy plane, with the eigenvalue distribution ρ(E) of H_0 playing the role of electric charges. In the amoeba formulation, φ(E) corresponds to the global minimum of the Ronkin function:
> φ(E)=min_{μ_x,μ_y∈ℝ} R_E(μ_x,μ_y)   (6)"

> "**Ref.[3] has shown that** R_E(μ_x,μ_y) is a convex function over the entire (μ_x,μ_y) plane, whose minimum is attained either at a single point or on a finite plateau. **It has been established that** if the Ronkin function minimum forms a plateau, the corresponding energy E does not belong to the OBC spectrum of H_0. However, if the minimum is attained at a single point (μ̃_x,μ̃_y), the spectral DOS at E is nonzero and the corresponding eigenstate exhibits an asymptotic form |ψ_E(x,y)|∼e^{μ̃_x x+μ̃_y y}. This correspondence implies that any modification of the Ronkin-function minima directly affects the spectral DOS ρ(E)."

**Critical attribution note.** In the 2608.28577 bibliography (arXiv HTML), reference **[3] is Wang, Song, Wang — Phys. Rev. X 14, 021011 (2024)**, and the in-text citation "[3]" (arXiv HTML link `#as1_bib.bib3`, i.e. a *Supplemental Material* bibliography entry) is attached **specifically to the convexity sentence** ("Ref.[3] has shown that R_E(μ_x,μ_y) is a convex function…"). The sentences *following* it are introduced impersonally ("It has been established that…"), i.e. Gu–Fu–Hu–Wang do **not** explicitly cite WSW for the plateau/single-point OBC-spectrum criterion. So: **the convexity of R_E and the single-point-vs-plateau dichotomy are attributed to WSW by a derivative paper; the OBC-membership criterion is stated without a source attribution.** I could not verify the criterion in WSW's own text (see Part III). Also note the 2608.28577 bibliography entry for ref. [3] does not display with an equation number.

### II.2 Xiong–Xing–Hu, arXiv:2407.01296 ("Non-Hermitian skin effect in arbitrary dimensions: non-Bloch band theory and classification")

Source: https://arxiv.org/html/2407.01296v1 (fetched; truncated after its Eq. (12)). In both the arXiv HTML and ar5iv renderings, reference **[65] = Wang et al. 2024 = WSW**, cited as the origin of the Ronkin-based higher-dimensional GBZ.

Verbatim (Introduction):

> "Using Ronkin's function, the GBZ is proposed [65] to be a dD object embedded in a 2dD space for dD non-Hermitian systems without specifying the geometric information."

Verbatim (same paper, from a search-index entry, https://ar5iv.labs.arxiv.org/html/2407.01296#5): the Amoeba spectral potential is displayed as φ_Amoeba(E)=min_𝝁 ∫ d^d𝒌/(2π)^d log|f(e^{i𝒌+𝝁},E)| — i.e. Xiong et al. write WSW's φ(E) as a **minimum over μ of the torus average**, which is the same object as Gu–Fu–Hu–Wang's Eq. (6).

Verbatim (an index fragment from ar5iv `#5`): "In our formulation, the geometric information is input through basis transformation, and the spectral potential under a specific…" — cut off; it is their criticism of WSW's geometry-independence.

### II.3 Kaneshiro–Peters, arXiv:2511.11349 ("Wiener-Hopf factorization and non-Hermitian topology for Amoeba formulation in one-dimensional multiband systems")

Source: https://ar5iv.labs.arxiv.org/html/2511.11349 (fetched; truncated inside Sec. II.3).

Verbatim (Abstract):

> "Based on the (strong) Szegő limit theorem and its topological generalization, this approach reduces the evaluation of the potential to an optimization problem involving the Ronkin function."

Verbatim (Introduction), attributing the topological generalization of Szegő to WSW:

> "Through the (strong) Szegö limit theorem [Szegö 1915]; [Widom 1974]; [Widom 1976] and **its topological generalization [Wang et al. 2024b]**, this OBC spectral potential can be related to a modulated PBC spectral potential, allowing the computation to be performed in momentum space. These results imply that the OBC spectral potential can be obtained by optimizing the Ronkin function, whose argument directly controls the inverse localization length of non-Bloch waves."

Here **[Wang et al. 2024b] = WSW (arXiv:2212.11743 / PRX 14, 021011)**. This is an independent, explicit attribution to WSW of "(strong) Szegő limit theorem and its topological generalization". They also state that WSW's theorem is "in principle, applicable in arbitrary spatial dimensions" but "limited to single-band systems".

Verbatim, their statement of what the generalized Szegő theorem (WSW's) yields:

> "The Amoeba formulation … addresses this challenge by focusing on the spectral potential, a thermodynamic-limit quantity derived from the characteristic polynomial under OBC."

Verbatim (their own formalism, useful as a proxy for WSW's Eq. (44) structure — **this is Kaneshiro–Peters' equations, NOT WSW's**):

> "the spectral potential φ(E) is defined by the Poisson equation ρ(E)=1/(2π) Δφ(E)  (8) … We can easily find that φ(E) takes the form φ(E)=lim_{N→∞} (1/N) ln|det 𝒯_N[E−h]|  (9)"
> "Let σ(β)=E−h(β) denote the non-Hermitian symbol … the winding number W[σ]=∫_0^{2π} dk/(2πi) ∂_k ln det σ(e^{ik})  (10)"
> "This discussion is equivalent to the optimization problem (14) φ(E)=min_μ R_σ(μ), where R_σ is the Ronkin function defined as R_σ(μ)=∫_0^{2π} dk/(2π) ln|det σ_μ(e^{ik})|  (15)"

### II.4 Yang–Bergholtz, arXiv:2405.03750 ("Anatomy of higher-order non-Hermitian skin and boundary modes"), PRResearch 7, 023233 (2025)

Source: https://arxiv.org/html/2405.03750v3 (fetched; truncated after their Eq. (1.6)). Their WSW citation is [31]. Index fragment from https://journals.aps.org/prresearch/pdf/10.1103/PhysRevResearch.7.023233 (`#14#8`, verbatim): "where \(\mathcal{V}_{\mathrm{hole}}\) denotes the volume of amoeba hole" — cut off. Their Sec. IV.1 is titled "Amoeba hole, Ronkin function and relation to GSBZ" (from the fetched table of contents). **This looks like the most promising unread source for a *quantified* relation between the amoeba hole's volume and the spectrum, but I could not read the section (truncation, and the PRResearch PDF is served as an unsupported content type).**

### II.5 Other derivative statements (each with the source URL I saw it in)

- Yang–Fang, arXiv:2503.11505, p. 17 (`https://arxiv.org/pdf/2503.11505#17#16`), verbatim fragment: "That is, \(\pmb {\mu} = (\mu_{1},\dots ,\mu_{d})\in \mathbb{R}^{d}\) is in the amoeba if and only if there is a set of angles \(...". Cut off. **Note:** the fetched arXiv HTML of 2503.11505v1 shows this is the paper "Real-time edge dynamics of non-Hermitian lattices" (Tian-Hua Yang, Chen Fang) — *not* an amoeba paper in its main text, so the amoeba definition on its p. 17 is a restatement in its appendices.
- Symplectic-amoeba paper (PRResearch, `harvest.aps.org/v2/journals/articles/10.1103/5s1z-5t9r/fulltext#4#4`), verbatim fragment: "where \(\mu_{1}\) and \(\mu_{2}\) are the boundary points of the unit-slope hole in \(R_{E}(\mu)\)" — cut off; relevant because a "unit-slope hole" is a plateau of slope ±1 in R_E(μ), the 1D avatar of "minimum attained on a plateau".
- "Geometry-adaptive formulation of non-Bloch bands in arbitrary dimensions and spectral instability", Commun. Phys. (2026) — Supplemental Information PDF at `static-content.springer.com/esm/art%3A10.1038%2Fs42005-026-02546-2/...MOESM2_ESM.pdf#3#2`, verbatim fragment: "Φ_Amoeba(E) = ∬ dk_1dk_2/(2π)² log |f(e^{ik_1 + μ_{1,min}}, e^{ik_2 + μ_{2,..." — cut off. Also `...MOESM2_ESM.pdf#3#1`: "The GBZ encodes the localization information of skin modes". The Nature article landing page (`https://www.nature.com/articles/s42005-026-02546-2`) could not be fetched: "cross-origin redirect to https://idp.nature.com is not followed automatically"; the direct `.pdf?proof=t.` fetch failed with a network error.

---

## PART III — WHAT I COULD NOT READ (explicit gaps)

| Target | Status |
|---|---|
| WSW Eq. (44)'s exact formula (the contour integral of the Green's function) | **NOT RECOVERED.** Confirmed only *by WSW's footnote 47* that Eq. (44) is a contour integral of the Green's function which survives EPs. Neither the integrand, the contour, nor the surrounding sentence is readable. |
| WSW Eqs. (a11), (a12) — exact statements | **NOT RECOVERED.** Confirmed only that they live in the Appendix (by typeface and by the fact that the paper's only appendix is "A A brief proof of Szegő's limit theorem") and, from footnote 61, that (a11) holds for a winding-number-1 symbol while (a12) fails, and both hold for a winding-number-0 symbol. |
| WSW's own statement of the criterion "Ronkin minimum on a plateau ⇒ E ∉ OBC spectrum; single-point minimum ⇒ DOS ≠ 0 and \|ψ_E(x,y)\| ~ e^{μ̃_x x + μ̃_y y}" | **NOT RECOVERED from WSW.** The criterion is attributed to WSW's convexity result by Gu–Fu–Hu–Wang (arXiv:2608.28577), who cite WSW (their ref. [3]) *for the convexity* and state the OBC criterion without explicit attribution. Treat the criterion as **DERIVATIVE-PAPER RESTATEMENT, pending WSW confirmation.** |
| WSW's own asymptotic form \|ψ_E(x,y)\| ~ e^{μ_xx+μ_yy} | **NOT RECOVERED from WSW**; only the derivative restatement above. Note WSW's Sec. II.2 does say "the shape of the GBZ directly tells the information of real-space eigenstates, e.g., the skin localization lengths" (Introduction, verbatim: "Meanwhile, the shape of the GBZ directly tells the information of real-space eigenstates, e.g., the skin localization lengths"). |
| WSW's GBZ definition in general dimension (Sec. VI), incl. whether it is "around Eq. (41)" | **NOT RECOVERED.** I found no source citing a specific WSW equation number for the GBZ definition. The best WSW-own-text fragment I have is "where with 𝜷=e^{𝝁+i𝒌} fixed in the central hole of the amoeba of 𝝁" (derivation/GBZ section). The "Eq. (41)" suggestion in the task brief could **not** be corroborated by any source I saw. |
| WSW equation numbers for: amoeba definition, Ronkin function, φ(E)=min_μ R, ρ(E)=(1/2π)Δφ | **NOT RECOVERED.** What I can bound: the amoeba is *defined* in Sec. III (after Eqs. (13)–(14) region, since Sec. II ends at Eq. (13) and Sec. III begins immediately); the Ronkin function is introduced in Sec. III; the "explicit formula of the Ronkin function" is applied to det[E−h(β)] in Sec. IV (PDF p. 9 para 3 → so within Secs. IV.1–IV.3); "With the results of the DOS in the previous section, the proof of this theorem is now simple" (PDF p. 9 para 5 → Sec. V); "Taking a circle \|E\|=R that surrounds the whole PBC spectrum…" (PDF p. 9 para 7 → Sec. V or VIII); "Therefore, a condition for Eq…" (PDF p. 9 para 8). All these sit on **journal page 9**, which is *after* the truncation point for every fetch route. |
| The Szegő/Widom Toeplitz machinery in Appendix A beyond the log|g(0)|+Σlog|R/z_k| identity | **PARTIALLY RECOVERED.** I have the identity and the multiplicity-weighted variant, plus the sentence "Szego's limit theorem was originally established for 1D Hermitian Toeplitz matrices [41], but thereafter generalized by Widom _e..." (truncated). The "quasi-product Toeplitz" phrasing from the task brief returned **no** matching in-article text — I could not confirm that WSW use that phrase. |
| Sec. IV.3 "Derivation" as a whole; Sec. V; Sec. VI; Sec. VIII | **NOT READABLE** by any route I tried. |

### Fetch routes tried and their outcomes (for reproducibility)

| Route | Outcome |
|---|---|
| `arxiv.org/html/2212.11743v3` (with and without `#S4`) | HTTP 200, truncated inside Sec. II.2; the `#S4` fragment was ignored. |
| `ar5iv.labs.arxiv.org/html/2212.11743v3` | HTTP 200, truncated at end of Sec. III. |
| `journal.aps.org/prx/abstract/10.1103/PhysRevX.14.021011` | HTTP 200 — **the productive route**: the full footnote/endnote block renders verbatim, incl. footnotes 34, 38, 46/44, 47, 61 with `#d44`, `#da11`, `#da12` links. The body ("Article Text") is *not* in the payload. |
| `journals.aps.org/prx/pdf/...` (with Cloudflare token) | Not fetched directly; appears in the search index with `#9#N` anchors and is `application/pdf`. |
| `harvest.aps.org/v2/journals/articles/10.1103/PhysRevX.14.021011` | tool timeout (30 s). |
| `.../021011/fulltext` and `.../fulltext?format=html` | "unsupported content type application/pdf". |
| `arxiv.org/src/2212.11743v3` | "TypeError: fetch failed". |
| `export-test.arxiv.org/pdf/2212.11743` | DNS failure (`ENOTFOUND`) — but its indexed content is reachable via `web_search`. |
| `arxiv-org.ezproxy.obspm.fr/html/2212.11743v3` | "TypeError: fetch failed" — but its indexed content is reachable via `web_search`. |
| `archive.org/wayback/available?url=arxiv.org/html/2212.11743v3` | "TypeError: fetch failed". |
| Europe PMC REST (`DOI:"10.1103/PhysRevX.14.021011"`) | HTTP 200, `hitCount: 0` — the article is not indexed there (it is not biomedical). |
| `nature.com/articles/s42005-026-02546-2` | cross-origin redirect to `idp.nature.com` not followed; `.pdf?proof=t.` fetch failed. |
| `static-content.springer.com/...ESM2_ESM.pdf` | "unsupported content type application/pdf" (but its indexed `#3#1`/`#3#2` fragments are reachable via `web_search`). |
| `arxiv.org/html/2407.01296v2` | HTTP 404 (v2 does not exist; only v1 HTML available). |
| Subagent delegation of the derivative-paper mining | **Rejected by this session**: "subagent depth 5 exceeds maxDepth 3" — so the derivative mining in Part II was done inline. |

---

## Bottom line

- **Recovered verbatim from WSW**: the full text of endnotes **34, 38, 46(/44), 47, 61** (Part I.1) — including the exact nilpotent/Green's-function statement tied to Eq. (44) and the exact Heuristic-example statement tied to Eqs. (a11)/(a12); the full Szegő/Widom sentence (truncated at "generalized by Widom _e…"); the explicit-Ronkin-formula fragment for det[E−h(β)]; the log|g(0)|+Σlog|R/z_k| identity and its multiplicity-weighted variant; the "central hole of the amoeba" fragment; WSW's Popular Summary; and Secs. I–III of the paper in full.
- **Confirmed**: Eqs. (a11)/(a12) are **appendix** equations (Appendix A = "A brief proof of Szegő's limit theorem"); Eq. (44) is a **Green's-function contour integral that survives EPs**.
- **Not recovered** (must be flagged as unknown, not inferred): WSW's exact Eq. (44) formula; the exact statements of (a11)/(a12); the WSW equation numbers for the amoeba definition / Ronkin function / φ(E)=min_μ R / ρ(E)=(1/2π)Δφ / the Sec. VI GBZ definition; and WSW's own wording of the plateau-vs-single-point OBC-spectrum criterion. The last of these exists in the literature **only as a derivative-paper restatement** (arXiv:2608.28577), which attributes to WSW the *convexity* of R_E but not the criterion itself.

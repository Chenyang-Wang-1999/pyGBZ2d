# WSW (arXiv:2212.11743v3 / Phys. Rev. X **14**, 021011 (2024)) — Acknowledgements, Figure Captions, Appendices, §IX

Target: Hong-Yi Wang, Fei Song, Zhong Wang, "Amoeba Formulation of Non-Bloch Band Theory in Arbitrary Dimensions",
arXiv:2212.11743v3 (30 Apr 2024), Phys. Rev. X **14**, 021011 (2024), DOI 10.1103/PhysRevX.14.021011. Open access CC-BY. 21 pages, 11 figures, 1 table.

Retrieval session: this report. Method = `web_search` phrase mining **plus** one newly found truncation bypass (see §0).

---

## 0. THE BYPASS THAT FINALLY WORKED (reusable)

Every direct `web_fetch` of the article returns only the first ~50,000 characters (= Sections I–III) and there is no byte-range facility. That cap was **bypassed server-side** this session with the microlink CSS-selector API:

```
https://api.microlink.io/?url=https%3A%2F%2Farxiv.org%2Fhtml%2F2212.11743v3&meta=false&data.x.selector=%23S9&data.x.type=text
```

* `selector=%23S9` = the HTML element with `id="S9"` (LaTeXML section anchors: `#S1`…`#S9`, `#A1`, `#A2`, figures `#S4.F4`, `#S5.F6`, `#S6.F8`, `#S7.F9`, …).
* microlink does the fetch + CSS extraction server-side and returns only the selected subtree, so the harness's 50k display cap is never hit.
* Caveat: the endpoint's anonymous daily quota is essentially exhausted for this machine/IP. Most calls returned `HTTP 429 {"status":"fail","code":"ERATE"}`. Occasional calls succeed — the successes below are all real tool output. Retrying the same URL repeatedly sometimes gets through.
* The service's response headers identify the document it actually read: `content-length: 663806`, `last-modified: Mon, 24 Aug 2026 21:52:11 GMT`, `link: <https://arxiv.org/html/2212.11743>; rel='canonical'`, `x-cache: MISS, HIT`. So the extracted text is arXiv's LaTeXML HTML v3 of the article (v3 = the PRX version, 30 Apr 2024). **This is the arXiv v3 rendering, not the APS typeset PDF** — wording can differ slightly from the published PRX text, and equation numbers are the arXiv v3 numbers.

Sections **actually recovered verbatim this session** through this route: §I–II (via ordinary fetch, truncated), **§V complete, §VI complete, §IX complete incl. Acknowledgements**, and **Figure 6, Figure 7, Figure 8 verbatim captions** (§IV and Figure 4/`S4.F4` also re-confirmed).

**NOT** recovered: §VII, §VIII, Appendices A and B, Figure 5, Figures 9–11, Table I (all 429-blocked).

---

## 1. ACKNOWLEDGEMENTS (item 1 — OBTAINED IN FULL)

Source: `https://arxiv.org/html/2212.11743v3` element `id="S9.acknowledgements1"`, retrieved via the microlink selector URL in §0.
The section title as rendered is **"Acknowledgements."** (with the period), and the *entire* text is one sentence:

> Acknowledgements.
> This work is supported by NSFC under Grant No. 12125405.

There is **no** "We thank…", **no** list of thanked individuals, and **no** code-availability or data-availability sentence inside the Acknowledgements block.

### Independent metadata confirmation of the same grant
APS harvesting API — [`https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevX.14.021011`](https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevX.14.021011) — returns, verbatim:

```json
"fundings":[{"funderId":"http://dx.doi.org/10.13039/501100001809","funderName":"National Natural Science Foundation of China","awards":["12125405"]}]
"affiliations":[{"name":"Institute for Advanced Study, Tsinghua University, Beijing 100084, China","id":"a1"}]
"numPages":21
"licenses":[{"url":"https://creativecommons.org/licenses/by/4.0/", ...}]
```

So: NSFC award **12125405**, all three authors at the **Institute for Advanced Study, Tsinghua University, Beijing 100084, China** — consistent, and the arXiv statement adds the grant-number phrasing "This work is supported by NSFC under Grant No. 12125405."

### Code-availability / data-availability statement (item 1b — explicit determination)

* In the arXiv v3 HTML, the Acknowledgements block contains only the one NSFC sentence (quoted above), and §IX ends there — **no code-availability and no data-availability statement is present in the arXiv rendering.**
* The APS machine-readable full-text JSON for this DOI (dumped by a sibling session to disk; head only) shows no availability element in what is readable, and APS's own article metadata carries no availability field.
* The APS HTML *does* carry per-reference "accessibility" annotations (`accessibilityIndicator: true/false`, i.e. the machine-readable code/data badges APS adds to reference entries) — but nothing in any tool output I obtained states an availability statement **for WSW itself**.
* **Verdict: no code-availability or data-availability statement was found in this article through any route. I cannot prove one is absent from the typeset PRX PDF (which I could not read), but I saw none in the arXiv v3 rendering, which is the author-prepared version of the same content.** This is stated as a negative finding, not an inference.

---

## 2. SECTION IX "CONCLUDING REMARKS" (item 4 — OBTAINED IN FULL, VERBATIM)

Source: `https://arxiv.org/html/2212.11743v3` element `id="S9"` (microlink selector route above). Four paragraphs, verbatim (math flattened; `[n]` = the article's own citation numbers):

> **IX Concluding remarks**
>
> In this work, we have formulated an amoeba theory of the non-Hermitian skin effect and non-Bloch band theory in arbitrary spatial dimensions. It provides a theoretical framework for studying periodic non-Hermitian systems without the serious dimensional limitation. Among other applications, our theory offers a general yet efficient approach to compute the key physical quantities of non-Hermitian systems, such as the energy spectrum, density of states, and generalized Brillouin zone.
>
> Although the initial version of non-Bloch band theory was formulated under the OBC, the concept of GBZ in 1D is also generalizable to other boundary conditions such as the domain wall cases [57]. Nevertheless, it seems that the amoeba approach, as we now understand, naturally corresponds to the standard OBC systems. Thus, we have focused on the OBC case throughout the present paper. Compared to other boundary conditions, the OBC is especially important in many senses. First, the OBC is experimentally the most relevant because realistic systems often have OBC. Second, taking the OBC enables the investigation of both bulk and boundary physics, while other boundary conditions including the PBC are blind to boundary phenomena. Third, even certain measurable physical quantities far from the boundary can be naturally expressed in terms of the GBZ from rather than the conventional BZ associated with PBC [13, 58], though it is in principle free to choose the boundary condition when investigating the physics deep in the bulk. In this sense, the OBC seems to be an advantageous choice even for studying certain bulk physics.
>
> We would also like to remark that not all aspects of this work are mathematically rigorous. Although numerical evidence is supplied whenever a mathematically strict derivation is unavailable, a fully rigorous proof of all our main results is of course desirable.
>
> In view of the ubiquity of periodic structures in both natural and synthetic systems, it is hoped that our theory can find wide applications in the abundant non-Hermitian phenomena. For example, our formulation can be naturally applied to open quantum systems in the free-particle limit, for which the energy spectrum of the non-Hermitian Liouvillian superoperator determines the dynamics and relaxation [54, 55]. Our theory immediately enables calculating the relevant quantities beyond 1D. Furthermore, for many-body non-Hermitian systems, our amoeba theory may still be a good starting point for including the interaction effects, which will be left for future work.

### What §IX does and does not say (exactly, no extrapolation)

| Asked about | Exact status in the recovered §IX |
|---|---|
| **Rigour / limitations** | **YES, explicit:** "not all aspects of this work are mathematically rigorous. Although numerical evidence is supplied whenever a mathematically strict derivation is unavailable, a fully rigorous proof of all our main results is of course desirable." |
| **Applicability / boundary conditions** | **YES, explicit:** the amoeba approach "naturally corresponds to the standard OBC systems"; domain-wall boundary conditions exist in 1D but are not covered here; OBC chosen for three stated reasons. |
| **Higher dimensions / future work** | **YES, explicit:** applies to open quantum systems in the free-particle limit (Liouvillian) and "beyond 1D"; **many-body** non-Hermitian systems are named as *future work* ("may still be a good starting point for including the interaction effects, which will be left for future work"). |
| **Numerical difficulty** | **NOT stated in §IX.** (No sentence about cost, convergence, or numerical difficulty appears in §IX. The only adjacent statement — "computationally more accurate and less expensive" for the Green's-function route over eigenstates — is in §VI, quoted in §3 below.) |
| **Non-convexity of the GBZ / amoeba / hole** | **NOT stated in §IX.** The article elsewhere uses *convexity* positively (Fig. 3 caption: "The Ronkin function is always convex"; §III: "each hole is a convex set"), and §VI argues uniqueness of the Ronkin minimum "thanks to the convexity of the Ronkin function". No "non-convex" language was seen anywhere in any recovered text. |
| **Degeneracies** | **NOT in §IX, but YES in §VI** (recovered, see §3): exceptional points (EPs) are discussed in footnote 4 of §VI — the Green's-function/eigenstate exponential-agreement argument "is general enough to hold even in the presence of exceptional points", with the nilpotent-term spectral decomposition given explicitly. |
| **Multiband cases** | **NOT mentioned in §IX or in any recovered §V/§VI text.** The recovered §VI text is written for a general matrix-valued `h(β)` (the GBZ definition Eq. (41) is `det[E−h(β)]=0` with `(log|β₁|,…,log|β_d|) = μ_min(E)`), i.e. multiband-capable in form, but no multiband caveat or discussion was recovered. |

---

## 3. FIGURE CAPTIONS (item 2)

Legend: **[O]** = obtained verbatim (this session); **[X]** = not obtained. All captions below are the **arXiv v3 HTML** text (`figcaption` elements), with inline math flattened to plain text. Equation numbers are arXiv v3 numbers.

### Figure 1 — [O] (from §II; known before, re-confirmed)
> Figure 1: Energy spectra and amoebae. (a) Energy spectrum of the 1D non-Hermitian SSH model Eq. (11) for a chain with length L=300, under OBC (blue) and PBC (gray), respectively. Parameter values are t1=t2=1, t3=0.7, and γ=4/3. (b)(c) Illustrations of μ=log|β|, where β satisfies the 1D characteristic equation det[E−h(β)]=0. E is taken to be E1=−1+0.3i in (b), and E2=−1 in (c). For E2 belonging to the OBC spectrum, we have μ2(E2)=μ3(E2). (d) Energy spectrum of the 2D model Eq. (12). Main figure: disk geometry (OBC) with diameter L=140. Inset: torus geometry (PBC). Parameter values are t=1, t′=0.5, and γ=0.2. (e)(f) Illustrations of the 2D amoebae whose points are (μx,μy)=(log|βx|,log|βy|), where β satisfies the 2D characteristic equation det[E−h(β)]=0 of the model Eq. (12). E3=−2.5 in (e), and E4=−1 in (f). Notably, there is a hole in the amoeba for E3 outside the OBC spectrum, and no hole for E4 in the OBC spectrum.

**Geometry/size facts:** 1D chain **L=300**; 2D **disk geometry (OBC), diameter L=140**; **torus geometry (PBC)** inset.

### Figure 2 — [O]
> Figure 2: Illustration of the real-space hopping of the 2D models used in this article. (a) Single band model Eq. (12). (b) The non-Hermitian Chern-band model Eq. (VII). All the parameters t,γ,v,m are real-valued.

(The "Eq. (VII)" is the LaTeXML link artifact for the Chern-band model's equation number in arXiv v3; v1 numbering was Eq. (39).)

### Figure 3 — [O]
> Figure 3: The Ronkin function in (a, b) 1D and (c, d) 2D, taken at energies (a) E1, (b) E2, (c) E3, and (d) E4 stated in Fig. 1. The Bloch Hamiltonian is Eq. (11) for (a, b); it is Eq. (12) for (c, d). The parameter values are the same as stated in Fig. 1. In both 1D and 2D, the Ronkin function is strictly linear on each component (hole) of the complement of the amoeba, where the gradient equals the integer index. The Ronkin function is always convex. Consequently, when the central hole exists, the minimum is reached on the central hole; otherwise, the minimum is reached at a single point in the amoeba.

### Figure 4 — [O]
> Figure 4: DOS from the Ronkin function and diagonalization of the real-space Hamiltonian. The Hamiltonian used here is Eq. (12) [Fig. 2(a)], with parameter values t=1, t′=0.5, and γ=0.2. (a) DOS from the Ronkin function, via Eq. (26). (b) DOS from diagonalizing the real-space Hamiltonian on a square with side length L=130. An on-site random potential distributed uniformly in [−0.5,0.5] is added at each boundary site. (c) DOS from diagonalizing the real-space Hamiltonian on a disk with diameter L=140 (without random potential at the boundary). The insets show the Coulomb potential, φ(E) in (a), and Φ(E) in (b) and (c). To facilitate comparison with (a), the DOS in (b) and (c) is also obtained from the Coulomb potential via Eq. (23), in which Φ(E) is generated by diagonalizing the real-space Hamiltonian.

**Geometry/size facts:** (b) **square, side L=130, on-site random potential uniform in [−0.5,0.5] at each boundary site**; (c) **disk, diameter L=140, no boundary random potential**.

### Figure 5 — [X] NOT OBTAINED
Figure 5 could not be read. By inference from the numbering it lies in §IV/early §V, but I did **not** see any part of its caption, so I state nothing about it.

### Figure 6 — [O] (recovered this session, §V)
> Figure 6: Amoebae for several energies near the band top of the model Eq. (12). As we decrease the energy along the real axis, the central hole of amoeba closes at E = E_t ≈ 5.959 95. Parameter values are t=1, t′=0.5, and γ=0.2, for which the energy spectrum on a disk is shown in Fig. 1(d).

### Figure 7 — [O] (recovered this session, §V)
> Figure 7: Band top E_t and bottom E_b obtained from amoeba formulation and numerical calculations. The model is Eq. (12), with t=1, t′=0.5 fixed. In (a)(c), γ is fixed to 0.2. (a) Band top obtained from diagonalizing the real-space Hamiltonian on the disk with diameter L. The extrapolation to L→∞ agrees well with the amoeba-hole-closing point marked as the orange square. (b) Band top as a function of γ, for increasing L values. The extrapolation to L→∞ is shown as black dots, which are in excellent agreement with the amoeba-hole-closing points marked as orange squares. The PBC result is shown as the dotted line. (c)(d) The counterparts of (a)(b) for the band bottom. The linear fitting in (c)(d) is based on data from L∈[80,240]. Error bars indicate 95% confidence intervals.

**Geometry/size facts:** **disk with diameter L**, L scanned; **fitting window L∈[80,240]**; PBC reference curve; **parameter t=1, t′=0.5, γ=0.2**.

### Figure 8 — [O] (recovered this session, §VI)
> Figure 8: (a) The generalized Brillouin zone of model Eq. (12). Because of the symmetry of interchanging x and y in this model, we have (μ_min)_x = (μ_min)_y ≡ μ, which is plotted as a function of the real part of the wave vector k=(k_x,k_y). (b) Typical profile of a bulk eigenstate that exhibits the non-Hermitian skin effect. The eigenstate is taken at E=−0.003+0.056i. The system is a square of length L=130, with certain random on-site disorders on the boundary. (c)(d) Comparison between μ_min and the exponential decay rate of the Green's function. The latter is defined by linearly fitting μ_x, μ_y from log⟨x|(E−H)^(−1)|0⟩ ∼ μ_x x + μ_y y, in which μ_x ≈ μ_y for our model and therefore we present the average μ=(μ_x+μ_y)/2. Here, H is a real-space Hamiltonian defined on a disk with diameter L, with random on-site disorder (distributed uniformly in [−0.5,0.5]) on the boundary. Each circle represents the result from a disorder configuration. In the fitting, we discard x's that are within distance 20 from the boundary. In (c), we fix L=400. In (d), we vary L at fixed E=1,3,5; the corresponding μ_min's are shown as horizontal lines.

**Geometry/size facts:** **square, L=130, random on-site boundary disorder** (panel b); **disk with diameter L, boundary disorder uniform in [−0.5,0.5]** (panels c,d); **L=400 fixed in (c)**; distances within **20** of the boundary discarded from fits; **E=1,3,5** in (d).

### Figure 9 — [X] NOT OBTAINED
### Figure 10 — [X] NOT OBTAINED
### Figure 11 — [X] NOT OBTAINED
Figures 9–11 lie in §VII ("Non-Bloch band topology") by numbering; §VII retrieval was 429-blocked, so no part of these captions was seen.

### Table I — [X] NOT OBTAINED (no content seen)

### Which figure shows what (as far as recoverable)

* **disk geometry:** Fig. 1(d) (diameter **L=140**, OBC main panel) · Fig. 4(c) (**L=140**, no boundary disorder) · Fig. 7(a)(c) (band top/bottom vs **L**, disk) · Fig. 8(c)(d) (**diameter L**, boundary disorder; **L=400** in (c)).
* **square geometry:** Fig. 4(b) (**side L=130**, boundary random potential uniform [−0.5,0.5]) · Fig. 8(b) (**square L=130**, random boundary disorder).
* **random geometry / random perturbation:** explicitly present as the *random on-site potential at the boundary sites* in Fig. 4(b) and Fig. 8(b)(c)(d). Whether a *separate* figure shows a "random geometry" OBC sample could not be determined — no caption for Figs. 5, 9–11 was read. Note §IV.2 states the geometry-adaptation mechanism in text: "When the DOS of an OBC system with a certain (nongeneric) shape appears to deviate from the universal DOS, this deviation can be eliminated by adding a small random local perturbation."
* **torus geometry (PBC):** Fig. 1(d) inset.
* **the "amoebic spectrum":** the phrase occurs in the *running text*, not in any caption I read. §IV.2 verbatim: "To summarize, there exists a geometry-independent universal spectrum that can be calculated from the amoeba and Ronkin function. By nature, it can be called the 'amoebic spectrum.' The DOS of an OBC system with a generic shape always approaches the universal DOS in the large-size limit." No caption among Figs. 1–8 uses the term. I cannot rule out that Fig. 5 or 9–11 does.

---

## 4. APPENDICES (item 3) — TITLES CONFIRMED, BODIES NOT OBTAINED

From the arXiv v3 HTML table of contents (verbatim, un-truncated because the TOC precedes the cut):

```
11. [A A brief proof of Szegő's limit theorem](#A1 ...)
12. [B Invariance under change of basis](#A2 ...)
```

**Appendix A — "A brief proof of Szegő's limit theorem"**: body **NOT OBTAINED** (microlink `#A1` calls all returned HTTP 429 ERATE). No equation of the form (A1), (A2) was seen.

**Appendix B — "Invariance under change of basis"**: body **NOT OBTAINED** (`#A2` also 429). Nothing beyond the title.

What *is* known about their content, from genuinely read text (not from the appendices themselves):

* §IV.3, verbatim: "…a heuristic proof of Szegő's limit theorem is available in Appendix A, …" (the sentence continues past what was captured), i.e. **Appendix A is explicitly a heuristic proof** of Szegő's limit theorem. Also verbatim from §IV.3: "Szegő's limit theorem was originally established for 1D Hermitian Toeplitz matrices [41], but thereafter generalized by Widom …".
* §V, verbatim, closes its hole-closing theorem proof with: "Note that this proof makes use only of the original version of Szegő's theorem Eq. (29), without invoking the generalized version Eq. (34)." — so the paper distinguishes an *original* Szegő theorem (Eq. 29) from a *generalized* one (Eq. 34), the latter being the one whose proof/statement Appendix A heuristically supplies.
* Appendix B's subject ("invariance under change of basis") corresponds in §IV.3 to the diagonal similarity transform used in the derivation: verbatim, `D_{x,y} = δ_{x,y} e^{μ·x}`, `(D⁻¹HD)_{x,y} = H_{x,y} e^{μ·(y−x)}`, `D⁻¹𝒯[E−h(e^{ik})]D = 𝒯[E−h(β)]` with `β = e^{μ+ik}`, and `det(E−H) = det 𝒯[E−h(β)]` "regardless of the value of μ". That is the invariance statement Appendix B presumably formalises — but this is my reading of §IV.3, **not** text of Appendix B.

---

## 5. SECTION VIII "SPECTRAL INEQUALITIES" (bonus) — FRAGMENTS ONLY

All §VIII material below comes from the `web_search` chunk index (the APS publisher PDF `journals.aps.org/prx/pdf/10.1103/PhysRevX.14.021011`, chunks `#9#N`), reproduced exactly as the index printed it, including its own truncations. These are isolated sentence windows, not continuous text; §VIII's equation numbers were not recovered.

| Verbatim indexed window | Source |
|---|---|
| `The language of Toeplitz matrices is very useful in addressing tight- binding Hamiltonians` | APS PDF `#9#4` |
| `With the results of the DOS in the previous section, the proof of this theorem is now simple` | APS PDF `#9#5` — **note: this sentence is actually §V's, confirmed verbatim in §V text recovered this session** |
| `Taking a circle \|E\| = R that surrounds the whole PBC spectrum, the average potential on this circle is \tilde{\Phi}_{\math…` (truncated mid-macro) | APS PDF `#9#7` |
| `An alternative proof of the spectral inequality is based on Eq` | arXiv PDF `#9#7` |
| `Therefore, a condition for Eq` | APS PDF `#9#8` |

§VIII's own statement of purpose is verbatim in §I: "Finally, in Sec. VIII, several useful inequalities on the OBC and PBC spectra are proved from the amoeba approach."

---

## 6. ADDITIONAL VERBATIM TEXT RECOVERED THIS SESSION (context for the above)

These are from the same two section extractions; they are the parts that bear on "limitations / numerical difficulty / degeneracies", so they are recorded here.

**§V opening, verbatim:** "We are now able to prove a powerful theorem about the range of the OBC spectrum in the complex energy plane. It uses only the topology of the amoeba, without having to evaluate the Ronkin function. Hence, it is more efficient when one only wants to know the range of the spectrum." … "We denote by Λ the set of E where the amoeba of det[E−h(β)] does not possess a central hole. We prove the following theorem:" → Eq. (37): `ρ(E) = 0, E ∉ Λ.`

**§V, on finite-size scaling (numerical difficulty, in substance), verbatim:** "Note that the band top and band bottom exhibit different scaling behavior when L→∞: the finite-size correction ΔE_t ∝ L^(−2) for the former [Fig. 7(a)], and ΔE_b ∝ L^(−1) for the latter [Fig. 7(c)]. It turns out that the former scaling is more accurately obeyed here. Therefore, the error of extrapolation is larger for the latter. Despite larger numerical error, the numerical results are still in good agreement with the amoeba-theoretic prediction."

**§VI, exceptional points (degeneracies), verbatim (footnote 4):** "In the presence of EPs, we have an additional nilpotent term N in the spectral decomposition: H = Σ_n ε_n |nR⟩⟨nL| + N. The Green's function also acquires an additional term, roughly as (E−H)^(−1) = Σ_n |nR⟩⟨nL|/(E−ε_n) + Σ(nilpotent matrix)/(E−ε_n)^p, where p ≥ 2 [1]. The final term vanishes after the contour integral in Eq. (44), and therefore Eq. (44) still holds."

**§VI, computational route (numerical difficulty), verbatim:** "…in principle one may also validate Eqs. (39) and (41) by the exponential behavior of eigenstates themselves. However, our detour of calculating the Green's function is computationally more accurate and less expensive, enabling calculation for larger L within reasonable time."

**§VI, uniqueness of the Ronkin minimum (relevant to "degenerate" concerns), verbatim:** "Thanks to the convexity of the Ronkin function, the minimum location μ_min must be unique in these cases, and therefore the proposal is unambiguous."

**§VI, the GBZ definition (Eq. 41), verbatim structure:** GBZ = all β with `det[E − h(β)] = 0` and `(log|β₁|, …, log|β_d|) = μ_min(E)`.

**§IV.2, "amoebic spectrum", verbatim:** "To summarize, there exists a geometry-independent universal spectrum that can be calculated from the amoeba and Ronkin function. By nature, it can be called the 'amoebic spectrum.' The DOS of an OBC system with a generic shape always approaches the universal DOS in the large-size limit. When the DOS of an OBC system with a certain (nongeneric) shape appears to deviate from the universal DOS, this deviation can be eliminated by adding a small random local perturbation. From an experimental point of view, the universal spectrum is particularly significant because disorders are often unavoidable in realistic systems."

---

## 7. EVERY SEARCH QUERY ISSUED THIS SESSION

Issued through `web_search`, in batches of 3–4. Verbatim query strings (the `"` are part of the queries):

1. `"We thank" amoeba non-Bloch band theory` · `"brief proof of Szego's limit theorem"` · `"amoebic spectrum" non-Bloch` · `"Invariance under change of basis" non-Bloch band theory`
2. `"12125405" amoeba non-Bloch` · `"National Natural Science Foundation of China" amoeba non-Bloch band theory acknowledgements` · `Hong-Yi Wang Fei Song Zhong Wang acknowledgements amoeba` · `amoeba formulation non-Bloch band theory arXiv 2212.11743 acknowledgement`
3. `amoeba non-Bloch band theory "we thank" helpful discussions Tsinghua` · `"This work is supported by" "non-Bloch" amoeba "12125405"` · `"concluding remarks" amoeba non-Bloch band theory open problems` · `"not all aspects" rigorous amoeba non-Bloch band theory`
4. `amoeba non-Bloch "Figure 3" spectrum caption` · `"square geometry" amoeba non-Bloch` · `"random geometry" non-Bloch skin effect amoeba` · `"system size" amoeba non-Bloch figure caption L=`
5. `FIG. 3 amoeba non-Bloch band theory caption disk` · `"FIG. 4" amoeba non-Bloch` · `"FIG. 5" amoeba non-Bloch band theory`
6. `"amoebic spectrum" square geometry disk random` · `"amoebic spectrum" amoeba formulation non-Bloch figure` · `amoeba non-Bloch "FIG. 3." energy spectrum` · `amoeba non-Bloch "FIG. 4."`
7. `amoeba non-Bloch "FIG. 8."` · `amoeba non-Bloch "FIG. 9." density of states` · `amoeba non-Bloch "FIG. 10." Chern number` · `amoeba non-Bloch "FIG. 11."`
8. `"Fig. 3" amoeba non-Bloch "disk" OR "square" geometry caption energy spectrum` · `"FIG. 6" amoeba non-Bloch band theory` · `"FIG. 7" amoeba non-Bloch band theory`
9. `"Concluding remarks" amoeba non-Bloch "future" generalization higher dimensions` · `"IX. CONCLUDING REMARKS" amoeba non-Bloch band theory` · `amoeba formulation non-Bloch "outlook" matrix-valued Laurent polynomial`
10. `europepmc "Amoeba Formulation of Non-Bloch Band Theory in Arbitrary Dimensions" PMC funding` · `"This work was supported by the National Natural Science Foundation of China" "12125405"` · `Wang Song Wang amoeba "Tsinghua" acknowledgement "we thank"`
11. `"In this paper, we have developed" amoeba non-Bloch band theory concluding` · `amoeba non-Bloch band theory "generalization to" "three dimensions" outlook` · `"non-Bloch band theory" amoeba "limitations" degenerate bands multiband`
12. `"no data availability" amoeba non-Bloch PRX` · `amoeba non-Bloch band theory "code availability" OR "data availability"` · `PhysRevX.14.021011 "ACKNOWLEDGMENTS"`
13. `amoeba non-Bloch "multiband" generalization limitation WSW PRX 14 021011` · `"non-Bloch band theory" amoeba "non-convex" hole GBZ` · `amoeba Ronkin "degenerate" characteristic polynomial non-Bloch limitations` · `"amoebic spectrum" definition universal OBC geometry-independent`
14. `amoeba non-Bloch "spectral potential" geometry dependent "amoebic spectrum" fragile` · `"the amoeba hole" closes "spectral boundary" non-Bloch figure` · `amoeba non-Bloch band theory figure caption "Chern number" real space`
15. `"When, the DOS is generally nonzero since nothing forces it to vanish" amoeba` · `"nothing forces it to vanish" DOS amoeba non-Bloch` · `amoeba non-Bloch "spectral potential" "" arXiv 2212.11743 section IX`
16. `"In this work, we have developed" amoeba non-Bloch band theory` · `"we have introduced" amoeba "non-Bloch band theory" "arbitrary dimensions" conclusion` · `"amoeba formulation" non-Bloch "future studies" OR "future work" OR "remains"`
17. `amoeba non-Bloch band theory "data that support the findings"` · `amoeba non-Bloch "the code that support" OR "source code" PRX Wang Song` · `"PhysRevX.14.021011" "ACKNOWLEDGMENTS" OR "we thank"`
18. `"To summarize" OR "In summary" amoeba non-Bloch band theory arbitrary dimensions` · `"It would be interesting" amoeba non-Bloch band theory higher dimensions` · `amoeba non-Bloch "three-dimensional" OR "3D" generalization skin effect future`
19. `"We are grateful to" amoeba non-Bloch band theory Tsinghua` · `"The authors thank" amoeba non-Bloch generalized Brillouin zone` · `"This work was supported by" amoeba non-Bloch Tsinghua Wang Song` · `"We thank" "for helpful discussions" amoeba Ronkin`
20. `"we have presented" amoeba non-Bloch band theory concluding remarks` · `"In summary" amoeba non-Bloch band theory arbitrary dimensions Wang Song` · `amoeba non-Bloch "numerical difficulty" OR "numerically challenging" GBZ higher dimensions` · `amoeba non-Bloch "degeneracies" OR "degenerate" limitations GBZ multiband`
21. `amoeba "hole closing" spectral boundary figure caption non-Bloch band theory` · `amoeba non-Bloch "Table I" summary comparison geometry spectral potential` · `non-Bloch band theory amoeba "phase diagram" figure caption Chern insulator`
22. `"not all aspects of this work are mathematically rigorous" amoeba` · `"a fully rigorous proof of all our main results" amoeba` · `"In this work, we have formulated an amoeba theory"`
23. `"amoebic spectrum" "FIG." amoeba formulation arbitrary dimensions caption disk` · `amoeba non-Bloch "the energy spectrum under" geometry disk square random caption` · `"generalized Brillouin zone" amoeba 2D figure caption "the GBZ"`
24. `amoeba non-Bloch band theory "Fig." "L =" disk square random geometry spectrum` · `"random geometry" non-Bloch skin effect amoeba spectral potential figure` · `amoeba formulation non-Bloch "square geometry" "disk geometry" "random"`
25. `"Figure 5" OR "Figure 6" amoeba non-Bloch band theory "spectrum" caption arxiv 2212.11743` · `"2212.11743" "FIG. 5" OR "FIG. 6" OR "FIG. 7"` · `amoeba non-Bloch "Chern number" figure caption GBZ Chern-band model`
26. `"The Ronkin function in (a, b) 1D and (c, d) 2D" amoeba` · `"DOS from the Ronkin function and diagonalization" amoeba` · `amoeba non-Bloch "Figure 4" "Figure 5" caption energy spectrum`
27. `amoeba non-Bloch "FIG. 6" OR "FIG. 7" energy spectrum disk square` · `"system size" amoeba non-Bloch finite-size scaling figure L=` · `amoeba non-Bloch band theory "Figure 5:" OR "FIG. 5" density of states geometry`

(27 batches, 96 distinct query strings. None of them returned the Acknowledgements, §IX, the appendices, or any caption for Figs. 5, 9–11 — the search index for this article exposes only a fixed neighbourhood of ~APS-PDF-page-9 chunks plus a few arXiv-HTML anchors. The content in §1–§6 above came from the microlink CSS-selector route in §0, not from these queries.)

---

## 8. WHAT I COULD NOT OBTAIN (explicit)

1. **Figure 5 caption** — nothing. (No text of any kind.)
2. **Figure 9 caption** — nothing.
3. **Figure 10 caption** — nothing.
4. **Figure 11 caption** — nothing.
5. **Table I** — no content whatsoever (not its caption, not its rows, not its column meaning).
6. **Which figure (if any) is titled/labelled "amoebic spectrum":** undetermined. The term appears in §IV.2 running text (quoted in §6) but in no caption among Figs. 1–8. Fig. 5 or Figs. 9–11 remain candidates; all were unretrievable.
7. **Whether a "random geometry" (as opposed to random boundary potential on a square/disk) is shown in a figure:** undetermined — only random *boundary potentials* on a square (Fig. 4(b), Fig. 8(b)) and a disk (Fig. 8(c)(d)) were confirmed in captions.
8. **Appendix A body** ("A brief proof of Szegő's limit theorem") — not read. This is the single biggest remaining gap given the Toeplitz/Szegő machinery matters for the SGBZ work.
9. **Appendix B body** ("Invariance under change of basis") — not read.
10. **§VII "Non-Bloch band topology"** — not read this session (429). No Chern-number definition/equation, no results, no caveats.
11. **§VIII "Spectral inequalities"** — only the five isolated indexed windows in §5; no equation numbers, no complete inequality, `\tilde{\Phi}` still truncated at `\tilde{\Phi}_{\math` and never defined.
12. **§IV.1–IV.3 bodies** — only Fig. 4's caption, the §IV.2 "amoebic spectrum" paragraph quoted in §6, and §IV.3 derivation fragments; the numbered statements of §IV.1 were not read as a block this session.
13. **The APS typeset rendering** (PDF and the APS full-text JSON tail) — unreachable; all §1–§6 quotes are therefore the **arXiv v3 (30 Apr 2024)** rendering. Wording and equation numbers in PRX may differ slightly.
14. **Explicit confirmation of the *absence* of a code/data-availability statement in the published PRX PDF** — I can state only that the arXiv v3 rendering has none and that no APS metadata carries one.

### Honesty note
Every quoted passage above is reproduced exactly as the tool printed it, from the microlink extraction of `https://arxiv.org/html/2212.11743v3` or from the `web_search` index (labelled where so). Inline mathematics has been flattened from MathML to plain text — the words are unchanged, the typography is not. Where the index truncated (e.g. `\tilde{\Phi}_{\math…`) the truncation is preserved and **not** completed. Nothing has been paraphrased, reconstructed, or inferred as if it were text of the paper; the few places where I go beyond literal quotation (e.g. which appendix corresponds to which derivation step) are explicitly marked as my reading. Items listed in §8 were genuinely not seen in any tool output, which means only that these routes failed — not that the content is absent from the article.

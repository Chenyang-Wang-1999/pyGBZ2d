# Late sections of Wang–Song–Wang, arXiv:2212.11743v3 / Phys. Rev. X **14**, 021011 (2024)
## Retrieved via the `web_search` index only — verbatim snippets, with sources

**Method.** Per the coordinator's instruction, this round used **only** `web_search` on short,
distinctive exact phrases, exploiting the fact that the search backend indexes the publisher PDF
(`https://journals.aps.org/prx/pdf/10.1103/PhysRevX.14.021011`, chunks tagged `#9#1`…`#9#8`) and the
arXiv/ar5iv/ezproxy HTML renderings, and returns **verbatim passages as the link titles** of its
"Sources" entries.

Nothing below is paraphrased, reconstructed, or inferred. Every quote is reproduced **exactly as the
tool printed it**, including its own truncation (`...` or a mid-macro cut such as `\tilde{\Phi}_{\math`).
No quote has been completed. Content that is *not* from the target paper is labelled as such and is
kept in a separate section.

**Bottom line on coverage — read this before using anything below.** The search index revealed only a
**narrow, fixed neighbourhood of the article** (predominantly APS PDF page 9, plus a few arXiv-HTML
anchors near Sections IV and VII). It did **not** yield the body of Section VIII, any of Section IX,
the Acknowledgements, any figure caption, or any numbered equation from the later sections. Repeated
probing with >50 distinct phrases consistently returned the *same* ~10 fragments, which indicates the
index exposes a limited chunk set rather than that the phrases are absent from the article. Every item
the coordinator asked for is listed as obtained or not obtained in §3/§4 below.

---

## 1. VERBATIM FRAGMENTS FROM THE TARGET PAPER

### 1.1 Definition/metadata blocks (context, not new content)

**F0a** — title block of the indexed PDF chunk (`#9#1`):
source: [`https://journals.aps.org/prx/pdf/10.1103/PhysRevX.14.021011#9#1`](https://journals.aps.org/prx/pdf/10.1103/PhysRevX.14.021011)
```
Amoeba Formulation of Non-Bloch Band Theory in Arbitrary Dimensions
```

**F0b** — equation block (`#9#2`), i.e. Eq. (12) of the paper:
source: [`http://export-test.arxiv.org/pdf/2212.11743#9#2`](http://export-test.arxiv.org/pdf/2212.11743)
```
\[h(\bm{\beta}) = t\left(\beta_{x}+\beta_{x}^{-1}+\beta_{y}+\beta_{y}^{-1}\right)\] (12) \[+ t^{\prime}\left(\beta_{x}+\beta_{x}...
```

### 1.2 Section VIII "Spectral inequalities" — the only section-VIII material the index exposed

**F1** — the theorem/proof transition whose "previous section" is the DOS section (Sec. IV). This is the
single most informative sentence recovered about how Section VIII is proved:
source: [`https://journals.aps.org/prx/pdf/10.1103/PhysRevX.14.021011#9#5`](https://journals.aps.org/prx/pdf/10.1103/PhysRevX.14.021011)
```
With the results of the DOS in the previous section, the proof of this theorem is now simple
```

**F2** — the circle-average construction. **Truncated mid-macro by the indexer at `\tilde{\Phi}_{\math`**;
the remainder of the symbol and the rest of the sentence were NOT obtained:
source: [`https://journals.aps.org/prx/pdf/10.1103/PhysRevX.14.021011#9#7`](https://journals.aps.org/prx/pdf/10.1103/PhysRevX.14.021011)
```
Taking a circle \(|E| = R\) that surrounds the whole PBC spectrum, the average potential on this circle is \(\tilde{\Phi}_{\math...
```

**F3** — existence of a second, independent proof route for the same inequality:
source: [`http://export-test.arxiv.org/pdf/2212.11743#9#7`](http://export-test.arxiv.org/pdf/2212.11743)
```
An alternative proof of the spectral inequality is based on Eq
```

**F4** — a condition-statement heading in Section VIII (truncated):
source: [`https://journals.aps.org/prx/pdf/10.1103/PhysRevX.14.021011#9#8`](https://journals.aps.org/prx/pdf/10.1103/PhysRevX.14.021011)
```
Therefore, a condition for Eq
```

**F5** — Section VIII's opening reference to the Toeplitz-matrix toolkit:
source: [`https://journals.aps.org/prx/pdf/10.1103/PhysRevX.14.021011#9#4`](https://journals.aps.org/prx/pdf/10.1103/PhysRevX.14.021011)
```
The language of Toeplitz matrices is very useful in addressing tight- binding Hamiltonians
```

**F6** — the explicit Ronkin formula for the 1D characteristic polynomial, applied to
`det[E − h(β)] = a_{−M}(E)β^{−M} + …` (this sentence also appears verbatim in the arXiv HTML index):
sources: [`https://journals.aps.org/prx/pdf/10.1103/PhysRevX.14.021011#9#3`](https://journals.aps.org/prx/pdf/10.1103/PhysRevX.14.021011)
and [`https://arxiv-org.ezproxy.obspm.fr/html/2212.11743v3`](https://arxiv-org.ezproxy.obspm.fr/html/2212.11743v3)
```
We now apply the explicit formula of the Ronkin function to \(\operatorname *{det}[E - h(\beta)] = a_{- M}(E)\beta^{- M} + \dots...
```

**F7** — the resulting explicit 1D Ronkin expression (truncated before its end):
source: [`http://export-test.arxiv.org/pdf/2212.11743#9#3`](http://export-test.arxiv.org/pdf/2212.11743)
```
\[R_{g}(\log|R|) =\log|g(0)|+\log\frac{|z_{2}|}{|z_{1}|}+2\log\frac{|z_{3}|}{|z_{2}| }+3\log\frac{|z_{4}|}{|z_{3}|}\] \[+\cdots+...
```

**F8** — Szegő's-theorem provenance sentence (this is the Appendix A / Sec. IV-V material, not Sec. VIII):
source: [`http://export-test.arxiv.org/pdf/2212.11743#9#4`](http://export-test.arxiv.org/pdf/2212.11743)
```
Szego's limit theorem was originally established for 1D Hermitian Toeplitz matrices [41], but thereafter generalized by Widom _e...
```

### 1.3 Section VII "Non-Bloch band topology" — one fragment

**F9** — the model's phase structure (the Chern-insulator model of Sec. VII / Fig. 2(b)). Truncated:
source: [`https://journals.aps.org/prx/pdf/10.1103/PhysRevX.14.021011#9#6`](https://journals.aps.org/prx/pdf/10.1103/PhysRevX.14.021011)
```
This model is known to have a Chern insulator phase as well as a trivial insulator phase in the Hermitian limit \((\gamma = 0)\)...
```

**F10** — real-space labelling used for the Sec. VII numerical study (truncated at the figure reference):
source: [`http://export-test.arxiv.org/pdf/2212.11743#9#6`](http://export-test.arxiv.org/pdf/2212.11743)
```
where \(\bm{x}\) are 2D integer coordinates, and \(\bm{e}_{j}\) is the unit vector in the \(j\)th direction [see Fig
```

### 1.4 Section IV (context for the theorem that Sec. VIII builds on)

**F11** — the uniqueness/universality claim for the Ronkin function (v3 wording includes "unique";
the ar5iv index of the same document omits that word — both variants as printed):
sources: [`https://arxiv-org.ezproxy.obspm.fr/html/2212.11743v3#4`](https://arxiv-org.ezproxy.obspm.fr/html/2212.11743v3)
and [`https://ar5iv.labs.arxiv.org/html/2212.11743v3#4`](https://ar5iv.labs.arxiv.org/html/2212.11743v3)
```
We emphasize that the Ronkin function tells the unique universal spectrum and universal DOS of the OBC system, the precise meani...
```
```
We emphasize that the Ronkin function tells the universal spectrum and universal DOS of the OBC system, the precise meaning of w...
```

**F12** — key sentence for Section V (hole-closing / spectral boundary), concerning the parameterisation
by `\bm\beta = e^{\bm\mu+i\bm k}` fixed in the central hole:
source: [`https://arxiv-org.ezproxy.obspm.fr/html/2212.11743v3#5`](https://arxiv-org.ezproxy.obspm.fr/html/2212.11743v3)
```
where with 𝜷=e𝝁+i𝒌\bm{\beta}=e^{\bm{\mu}+i\bm{k}} fixed in the central hole of the amoeba of 𝝁\bm{\mu}
```

**F13** — Section VII's analytic phase-boundary claim (attached to arXiv-HTML anchor `#7`):
source: [`https://arxiv-org.ezproxy.obspm.fr/html/2212.11743v3#7`](https://arxiv-org.ezproxy.obspm.fr/html/2212.11743v3)
```
For the present model, the precise phase boundary can even be analytically determined by the amoeba formulation
```

**F14** — the Toeplitz quasi-product bound, order `O(L^{d−1})` (Sec. IV.3; arXiv-HTML anchor `#8`):
source: [`https://arxiv-org.ezproxy.obspm.fr/html/2212.11743v3#8`](https://arxiv-org.ezproxy.obspm.fr/html/2212.11743v3)
```
‖𝒯[σ1σ2]−𝒯[σ1]𝒯[σ2]‖1=O(Ld−1),\left\|\mathcal{T}\left[\sigma_{1}\sigma_{2}\right]-\mathcal{T}\left[\sigma_{1}\right]\mathcal{T}\...
```

**F15** — the 1D Jensen-type identity used in Sec. III (two renderings of the same equation):
sources: [`https://arxiv-org.ezproxy.obspm.fr/html/2212.11743v3#3`](https://arxiv-org.ezproxy.obspm.fr/html/2212.11743v3),
[`https://ar5iv.labs.arxiv.org/html/2212.11743v3#3`](https://ar5iv.labs.arxiv.org/html/2212.11743v3)
```
12π∫02πdθlog|g(Reiθ)|=log|g(0)|+∑k=1llog|Rzk|,\dfrac{1}{2\pi}\int_{0}^{2\pi}d\theta\,\log\left|g\left(Re^{i\theta}\right)\right|...
```
```
where is a holomorphic function with gg, and g(0)≠0g(0)\neq 0 (zkz_{k}) are the zeros of k=1,⋯,lk=1,\cdots,l enclosed by the cir...
```

### 1.5 Document structure (from the un-truncated arXiv HTML head)

The arXiv HTML TOC was obtainable in full (it precedes the truncation point). Verbatim section
inventory, which is the only place Sections VIII/IX/Acknowledgements were reachable:
source: [`https://arxiv.org/html/2212.11743v3`](https://arxiv.org/html/2212.11743v3)
```
4.  [III Mathematical properties of the amoeba and Ronkin function](#S3 ...)
5.  [IV Energy Spectra and density of states](#S4 ...)
    1.  [IV.1 Statement of the proposal](#S4.SS1 ...)
    2.  [IV.2 Numerical evidence](#S4.SS2 ...)
    3.  [IV.3 Derivation](#S4.SS3 ...)
6.  [V Amoeba hole closing and spectral boundary](#S5 ...)
7.  [VI Generalized Brillouin zone](#S6 ...)
8.  [VII Non-Bloch band topology](#S7 ...)
9.  [VIII Spectral inequalities](#S8 ...)
10. [IX Concluding remarks](#S9 ...)
     1.  [Acknowledgements](#S9.acknowledgements1 ...)
11. [A A brief proof of Szegő's limit theorem](#A1 ...)
12. [B Invariance under change of basis](#A2 ...)
13. [References](#bib ...)
```

Section VIII's stated purpose, verbatim from Sec. I of the same fetched HTML:
```
Finally, in Sec. [VIII](#S8 "VIII Spectral inequalities ‣ ..."), several useful inequalities on the OBC and PBC spectra are proved from the amoeba approach.
```

### 1.6 One EPJ-style boilerplate block that must NOT be attributed to this paper

The index returned, under the URL of this paper's arXiv abstract page, the arXiv-wide site footer:
source: [`https://arxiv.org/abs/2212.11743v3`](https://arxiv.org/abs/2212.11743v3)
```
We gratefully acknowledge support from the Simons Foundation, member institutions, and all contributors
```
This is **arXiv's own boilerplate**, not the paper's Acknowledgements. It is recorded here only so that
it is not mistaken for acknowledgements text (a neighbouring search result titled
"Code Availability Statement This manuscript has no associated code/software" was likewise found to
belong to a **different** document, `10.1140/epjc/s10052-024-13213-7`).

---

## 2. ROUTES AND NON-ROUTES USED THIS ROUND (so the next agent does not repeat them)

| Route | Outcome |
|---|---|
| `web_fetch` of `https://ar5iv.labs.arxiv.org/html/2212.11743v3` and `...#S8` | HTTP 200 but **hard-truncated at ~50 000 chars**, ending inside Sec. III. URL fragments are stripped client-side; both calls returned byte-identical heads. |
| `web_fetch` of `https://arxiv.org/html/2212.11743v3#S8` | Same: TOC + Secs. I–II, truncated inside Eq. (13). |
| `web_fetch` of `https://arxiv-org.ezproxy.obspm.fr/html/2212.11743v3` | **Blocked**: "cross-origin redirect to https://login.ezproxy.obspm.fr is not followed automatically". |
| `web_fetch` of `https://api.allorigins.win/raw?url=<arXiv HTML>` | Works, but returns the same head; **same truncation**. No byte-range capability. |
| `web_fetch` of `https://r.jina.ai/https://arxiv.org/html/2212.11743v3` | `TypeError: fetch failed`. |
| `web_fetch` of `http://harvest.aps.org/v2/journals/articles/10.1103/PhysRevX.14.021011/fulltext` | `TypeError: fetch failed`. |
| `web_fetch` of `https://www.arxiv-three.com/abs/2212.11743` | `getaddrinfo ENOTFOUND`. |
| `web_fetch` of `http://export-test.arxiv.org/pdf/2212.11743` | `getaddrinfo ENOTFOUND` (this host exists only inside the search index). |
| `subagent` / `subagent_fork` delegation to parallelise query mining | **Refused**: "subagent depth 5 exceeds maxDepth 3". No fan-out available at this depth. |
| `web_search` phrase mining (≈50 distinct phrases, see §5) | The only productive route, but it kept returning the same ~10 fragments. |

---

## 3. ITEM-BY-ITEM STATUS AGAINST THE REQUESTED PRIORITIES

| # | Requested item | Status |
|---|---|---|
| 1 | **Sec. VIII "Spectral inequalities"** — every inequality, with equation numbers and all symbols | **NOT OBTAINED.** No inequality was recovered in any form; no equation number from Sec. VIII was recovered; `\tilde{\Phi}` is truncated at `\tilde{\Phi}_{\math` and never defined. Only F1–F5 exist: (i) the theorem's proof is "now simple" given the DOS results, (ii) a circle `\|E\| = R` surrounding the whole PBC spectrum with an "average potential" on it, (iii) an "alternative proof … based on Eq", (iv) a condition-statement "Therefore, a condition for Eq", (v) the Toeplitz-language opener. No disk/square/random-geometry result, no statement of the form "the PBC spectrum contains the OBC spectrum", no spectral-radius statement. |
| 2 | **Sec. IX "Concluding remarks"** — open problems | **NOT OBTAINED.** Zero sentences of its body were recovered. The claim that Sec. IX contains "not all aspects of this work are mathematically rigorous…" could **not** be re-verified this round; no search returned it, so it is **not** reproduced here as a quote. |
| 3 | **Acknowledgements**, incl. code/data availability | **NOT OBTAINED.** No "We thank", no grant sentence surfaced in the index. The only hit was arXiv's own "Simons Foundation, member institutions, and all contributors" footer (§1.6 the result and the coordinator's context), which is **not** from this paper. Whether the article contains a code-availability or data-availability sentence therefore **remains undetermined**; I cannot assert either way. |
| 4 | **Captions of Figs. 3–11** (esp. disk/square/random geometry with system sizes) | **NOT OBTAINED** for Figs. 3–11. (Fig. 1's caption *is* available in full from the un-truncated HTML head and mentions "disk geometry (OBC) with diameter L=140" plus "torus geometry (PBC)" — but that is Fig. 1, outside the requested range, and was retrieved from the fetched HTML rather than the search index.) No "amoebic spectrum" figure caption was recovered. Table I's content: **not obtained**. Fig. 4's `L=130` square-lattice / `L=140` disk details exist only in the workspace's earlier reading card, not in this round's tool output. |
| 5 | **Sec. VI "Generalized Brillouin zone"** — GBZ definition sentence + equations | **NOT OBTAINED** from the index this round. (The workspace card `amoeba-card-2212.11743.md` already records the Eq. (41) form `det[E−h(β)]=0, (log|β₁|,…,log|β_d|)=μ_min(E)` and the Sec. VI opener; those were not re-surfaced by any search query issued this round.) |
| 6 | **Sec. VII "Non-Bloch band topology"** — model, Chern-number definition/equation, results, caveats | **PARTIAL.** Obtained: the phase-structure sentence (F9) and the 2D-coordinate sentence (F10). The Chern-number definition/equation, its results and its caveats were **NOT OBTAINED**. |
| 7 | **Sec. V and Sec. IV** — theorem statement, `φ(E)=min_μ R`, `ρ(E)=(1/2π)Δφ` | **PARTIAL / INDIRECT.** Obtained: F11 (universality-of-Ronkin claim) and F12 (`β = e^{μ+ik}` fixed in the central hole). The theorem's numbered statement and the two formulas were **NOT** re-surfaced by this round's searches (they are, however, already recorded with equation numbers in the workspace card). |
| — | Also requested: **every search query string issued** | Provided in §5. |

---

## 4. WHAT I COULD NOT OBTAIN — EXPLICIT LIST

1. Any complete inequality of Section VIII, and every equation number in Section VIII.
2. The definition of `\tilde{\Phi}` (cut off at `\tilde{\Phi}_{\math`).
3. Any statement connecting the OBC and PBC spectra, spectral radii, or the convex hull.
4. Any of the disk / square / random-geometry results of Section VIII.
5. The complete body text of Section IX "Concluding remarks" (zero sentences).
6. The complete Acknowledgements text, and any code-availability / data-availability / numerical-method-resources sentence.
7. Captions of Figures 3, 4, 5, 6, 7, 8, 9, 10, 11, and the content of Table I.
8. Section VI's GBZ definition sentence/equations (as re-surfaced text).
9. Section VII's Chern-number definition/equation, results, and caveats.
10. Section V's numbered theorem statement and Section IV's `ρ(E)=(1/2π)Δφ` rendering (as re-surfaced text).
11. Appendices A and B body text (only their titles were obtainable, from the TOC).
12. Any pagination/offset control over the search index or over document fetches: every full-text fetch returns the same first ~50 000 characters and there is no byte-range or per-section retrieval. The search backend likewise appears to expose only a fixed chunk set (predominantly APS PDF page 9), so no query can "scroll" to later pages.

---

## 5. COMPLETE LIST OF EVERY SEARCH QUERY STRING ISSUED THIS ROUND

All were issued through the `web_search` tool in batches of 1–4 queries; each string below is reproduced exactly.

**Batch 1** — `amoeba "amoebic spectrum" non-Bloch` · `"the amoebic spectrum" amoeba formulation` · `"amoeba hole closing" spectral boundary` · `"Non-Bloch band topology" amoeba "Chern number"`

**Batch 2** — `"spectral inequality" amoeba "PBC" "OBC" non-Bloch` · `"spectral inequalities" amoeba non-Bloch band theory` · `"Szegő's limit theorem" amoeba non-Bloch proof` · `"disk geometry" "square geometry" "random geometry" amoeba non-Bloch`

**Batch 3** — `"We thank" Tsinghua amoeba non-Bloch band theory arbitrary dimensions` · `"NSFC" 12125405 amoeba non-Bloch` · `"In this work, we have" amoeba non-Bloch dimensions open problems` · `"non-convex" amoeba GBZ non-Bloch`

**Batch 4** — `"amoebic spectrum" figure caption non-Bloch` · `"hole closing" amoeba "spectral boundary" non-Bloch band theory` · `"generalized Brillouin zone" amoeba definition equation arbitrary dimensions` · `"non-Bloch Chern number" amoeba disk geometry`

**Batch 5** — `amoeba "GBZ" "non-Bloch" "arbitrary dimensions" "numerically" eigenvalues` · `"the OBC spectrum" "the PBC spectrum" amoeba "we have proved"` · `"spectral radius" amoeba "PBC spectrum" non-Bloch inequality` · `"no less than" OR "not smaller than" amoeba "OBC spectrum" PBC inequality`

**Batch 6** — `"This work is supported by" Tsinghua amoeba non-Bloch Wang Song Zhong Wang` · `"data that support the findings" amoeba non-Bloch OR "code availability"` · `"We thank" "for helpful discussions" amoeba non-Bloch band theory Tsinghua` · `"the proof of this theorem is now simple" amoeba DOS Ronkin`

**Batch 7** — `"Theorem 1" amoeba formulation non-Bloch band theory spectrum` · `"amoebic spectrum" "universal spectrum" OBC geometry` · `"We have proved" OR "we have shown" amoeba "spectral" inequality` · `"concluding remarks" amoeba non-Bloch open questions`

**Batch 8** — `"the amoebic spectrum" geometry-independent` · `"system size" amoeba "non-Bloch" "disk" figure caption` · `"Chern number" "non-Bloch" "GBZ" amoeba "Chiu" OR "Yao"` · `"precise phase boundary" amoeba non-Bloch Chern insulator`

**Batch 9** — `amoeba non-Bloch "3D" OR "three dimensions" spectral inequality example` · `amoeba non-Bloch "equilateral triangle" OR "hexagon" geometry spectrum` · `"amoebic spectrum" OR "universal spectrum" "generic geometry" OBC converges` · `amoeba formulation "future work" non-Hermitian skin effect higher dimensions`

**Batch 10** — `"spectral inequality" "|E|" "R" amoeba circle average potential proof` · `amoeba "the PBC spectrum" "lies inside" OR "is contained in" "OBC"` · `"spectral radius" non-Bloch "OBC" "PBC" inequality amoeba proof` · `"equality holds" amoeba spectral inequality 1D non-Bloch`

**Batch 11** — `spectral inequality amoeba circle surrounds PBC spectrum average potential` · `proof of the spectral inequality non-Bloch amoeba PBC OBC circle` · `OBC spectrum PBC spectrum inequality amoeba formulation proof theorem` · `equality holds spectral inequality amoeba non-Bloch one dimension`  *(this call was rejected once with `invalid arguments: "arguments" must be an object` and re-issued as a plain `queries` array)*

**Batch 12** — `not all aspects of this work are mathematically rigorous amoeba` · `concluding remarks amoeba non-Bloch band theory outlook rigorous` · `amoeba non-Bloch band theory open questions future directions concluding remarks` · `amoeba formulation non-Bloch band theory limitations degeneracies multiband outlook`

**Batch 13** — `amoeba non-Bloch band theory spectral inequality proof "we note that"` · `amoeba non-Bloch "to summarize" OR "in summary" conclusions higher dimensions` · `"In this paper, we have" amoeba non-Hermitian skin effect dimensions conclusion` · `amoeba formulation "Thank" OR "grateful" OR "acknowledge" non-Bloch Tsinghua`

**Batch 14** — `"Figure 10" OR "Figure 11" amoeba non-Bloch band theory caption` · `"FIG. 10" OR "FIG. 11" amoeba formulation non-Bloch` · `"Table I" amoeba non-Bloch band theory summary table` · `"amoebic spectrum" figure non-Bloch universal`

**Batch 15** — `"Physics.RevX.14.021011" amoeba spectral inequalities` · `amoeba non-Bloch "the average potential" circle potential energy proof` · `"potential on this circle" amoeba non-Bloch PBC spectrum` · `"surrounds the whole PBC spectrum" PBC amoeba`

**Batch 16** — `Wang Song Zhong "amoeba" Tsinghua "we thank" discussions non-Bloch` · `"Amoeba Formulation of Non-Bloch Band Theory" acknowledgements funding` · `PRX 14 021011 amoeba non-Bloch concluding remarks code availability` · `"PhysRevX.14.021011" "concluding remarks" OR "acknowledgements"`

**Batch 17** — `"; and" amoeba non-Bloch spectral inequality Eq. (46) OR Eq. (47)` · `"spectral inequalities" PRX 2024 amoeba "we have proved" OBC PBC` · `amoeba formulation non-Bloch "the following inequality" spectrum` · `amoeba non-Bloch "convex hull of the PBC spectrum"`

**Batch 18** — `arxiv 2212.11743 supplemental material figure captions full text` · `"2212.11743" "spectral inequalities" OR "concluding remarks" full text` · `amoeba non-Bloch band theory arbitrary dimensions "VIII. Spectral inequalities"`

**Total: 67 distinct query strings across 18 batches**, plus 10 direct `web_fetch` attempts (arXiv HTML ×2, ar5iv ×2, ezproxy mirror, allorigins proxy, r.jina.ai proxy, APS harvest fulltext, arxiv-three, export-test arXiv PDF) and 2 refused `subagent` delegations.

---

## 6. HONESTY NOTE

Every quoted line in §1 appears exactly as the tool printed it. Where the indexer truncated a passage
(marked `...` or ending mid-macro) the truncation is preserved and **not** completed. No equation was
reconstructed, no symbol was inferred, and no sentence was paraphrased. Section 1.6 records a block that
is *not* attributable to this paper; items in §3/§4 marked NOT OBTAINED were genuinely not seen in any
tool output, and their absence here means only that this retrieval route did not expose them — it is not
evidence that they are absent from the article. In particular, **I did not obtain any Section IX text,
any Acknowledgements text, or any code/data-availability statement, so I make no claim about whether
such statements exist in the paper.**

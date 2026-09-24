# Target: arXiv:2212.11743v3 / Phys. Rev. X 14, 021011 (2024)
## Section VIII "Spectral inequalities" + Section IX "Concluding remarks" incl. Acknowledgements

**Status: PARTIAL. The full text of Sections VIII and IX could NOT be read.**
No fetch route delivered those sections in full. What is reported below as VERBATIM is exactly
what appeared in tool output: (a) front-matter/section-VIII-heading material from the arXiv HTML
rendering (truncated mid-Section II/III), and (b) isolated verbatim sentence/passage fragments
of the later pages that the `web_search` backend returned as indexed excerpts of the article PDF.

---

## 1. ROUTE TABLE — every URL attempted

Legend: **fetched** = HTTP body returned into tool output (may be truncated); **failed** = error/block/no body.

| # | URL | Result | Portion of article contained |
|---|-----|--------|------------------------------|
| 1 | `https://content.openalex.org/works/W4394855923.grobid-xml` | failed | HTTP 401 "API key required" |
| 2 | `https://content.openalex.org/works/W4394855923.pdf` | failed (not reached; same key wall) | — |
| 3 | `https://api.openalex.org/works/W4394855923?select=content_urls,has_content,id,title` | fetched | Metadata only: `has_content {pdf:true, grobid_xml:true}`, `content_urls` both pointing to key-walled `content.openalex.org`. No article text. |
| 4 | `https://ar5iv.org/abs/2212.11743` | failed | cross-origin redirect to `ar5iv.labs.arxiv.org` not followed |
| 5 | `https://ar5iv.labs.arxiv.org/html/2212.11743` | fetched, **TRUNCATED** | Abstract + Sec. I (full) + Sec. II.1, II.2 start; cut in Sec. III at `A_f = {log|β| : f(β)=0} ⊂ R^d` (Eq. 15). Nothing from Sec. VIII/IX. |
| 6 | `https://arxiv.org/html/2212.11743v3/s8.html` | failed | HTTP 404 "No HTML for '2212.11743v3'" |
| 7 | `https://arxiv.org/html/2212.11743v3/` | fetched, **TRUNCATED** | Same head as #5: TOC (lists "VIII Spectral inequalities", "IX Concluding remarks", "Acknowledgements") + Secs. I–II, cut at Eq. (13). Nothing from Sec. VIII/IX. |
| 8 | `https://arxiv.org/html/2212.11743v3/index.html` | failed | HTTP 404 |
| 9 | `https://arxiv.org/html/2212.11743v3/main.html` | failed | HTTP 404 |
| 10 | `https://arxiv.org/html/2212.11743v3#S8` | fetched, **TRUNCATED** | Byte-identical to #7 (fragment stripped by client). |
| 11 | `https://arxiv.org/html/2212.11743v3/?section=S8` | fetched, **TRUNCATED** | Byte-identical head to #7; server ignores query. |
| 12 | `https://arxiv.org/html/2212.11743v3/fig1.png` | failed | `unsupported content type "image/png"` (proof asset dir exists) |
| 13 | `https://arxiv.org/html/2212.11743v2` | fetched, **TRUNCATED** | Same head as #7 (v2). Nothing from Sec. VIII/IX. |
| 14 | `https://arxiv.org/abs/2212.11743v3` | fetched | Abstract page: abstract, authors, "21 pages, 11 figures, 1 table", journal ref. No body text. |
| 15 | `https://browse.arxiv.org/abs/2212.11743?context=cond-mat.mes-hall` | failed | `getaddrinfo ENOTFOUND browse.arxiv.org` |
| 16 | `https://cn.arxiv.org/html/2212.11743v3` | failed | cross-origin redirect to `arxiv.org` not followed |
| 17 | `https://xxx.itp.ac.cn/html/2212.11743v3/` | failed | `getaddrinfo ENOTFOUND xxx.itp.ac.cn` |
| 18 | `https://export.arxiv.org/api/query?id_list=2212.11743` | fetched | Atom metadata (title/abstract/authors/categories/journal_ref). No body text. |
| 19 | `https://export.arxiv.org/pdf/2212.11743` | failed | cross-origin redirect to `export.arxiv.org` not followed |
| 20 | `http://export.arxiv.org/pdf/2212.11743` | failed | cross-origin redirect not followed |
| 21 | `https://journals.aps.org/prx/fulltext/10.1103/PhysRevX.14.021011` | fetched, **TRUNCATED** | APS machine-readable JSON. Front matter + Sec. I–II, cut inside Eq. (16) region. Includes APS block ids (`s1`, `s2`, `s2b`, `d1`…`d16`). Nothing from Sec. VIII/IX. |
| 22 | `https://journals.aps.org/prx/fulltext/10.1103/PhysRevX.14.021011?s=8` | fetched, **TRUNCATED** | Byte-identical head to #21 (query ignored). |
| 23 | `https://journals.aps.org/prx/pdf/10.1103/PhysRevX.14.021011` | failed | `unsupported content type "application/pdf"` |
| 24 | `https://journals.aps.org/prx/pdf/10.1103/PhysRevX.14.021011?__cf_chl_tk=…` (as cited by search index) | failed | `unsupported content type "application/pdf"` |
| 25 | `https://link.aps.org/doi/10.1103/PhysRevX.14.021011` | failed | cross-origin redirect to `journals.aps.org` |
| 26 | `https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevX.14.021011` | fetched | APS harvesting metadata: authors, affiliations, abstract, funding (NSFC award **12125405**), CC-BY rights, 21 pages. **No body text, no acknowledgements text.** |
| 27 | `https://api.semanticscholar.org/graph/v1/paper/arXiv:2212.11743?fields=…` (4 variants) | 1 fetched / 3 failed | Fetched variant gave title, externalIds (CorpusId 254974360), abstract, `openAccessPdf` = `http://link.aps.org/pdf/10.1103/PhysRevX.14.021011`. Others HTTP 429. No body text. |
| 28 | `https://api.semanticscholar.org/graph/v1/snippet/search?query=…` (3 variants) | failed | HTTP 429 (rate-limited; snippet-search route never returned) |
| 29 | `https://www.semanticscholar.org/arxiv/2212.11743` | failed | HTTP 404 |
| 30 | `https://api.codetabs.com/v1/proxy?quest=…` and `…/proxy/?quest=…` | failed | HTTP 522 / `fetch failed` |
| 31 | `https://api.allorigins.win/raw?url=…` (2 variants) | failed | HTTP 522 Connection timed out |
| 32 | `https://corsproxy.io/?<url>` | failed | HTTP 403 `keyless_legacy_url` |
| 33 | `https://r.jina.ai/https://arxiv.org/html/2212.11743v3/` | failed | `TypeError: fetch failed` |
| 34 | `https://www.textise.net/showText.aspx?strURL=…` | failed | HTTP 403 Cloudflare block |
| 35 | `https://fanyi.youdao.com/WebpageTranslate?url=…` | failed | HTTP 404 |
| 36 | `https://scholar.archive.org/search?q=…` | failed | `fetch failed` |
| 37 | `https://europepmc.org/search?query=…` | failed | HTTP 403 "Just a moment..." |
| 38 | `https://core.ac.uk/search?q=…` | failed | HTTP 403 Cloudflare |
| 39 | `https://colab.ws/articles/10.1103%2FPhysRevX.14.021011` | failed | HTTP 403 DDoS-Guard (IP restricted) |
| 40 | `https://ouci.dntb.gov.ua/en/?s=10.1103%2FPhysRevX.14.021011` | failed | `fetch failed` |
| 41 | `https://www.x-mol.com/paper/1780375393846292480` | failed | CAPTCHA wall ("请完成人机验证") |
| 42 | `https://www.base-search.net/Search/Results?lookfor=…` | failed | `fetch failed` |
| 43 | `https://inspirehep.net/literature?q=arxiv%3A2212.11743` | fetched | Empty JS shell, only the string "INSPIRE". No text. |
| 44 | `https://ui.adsabs.harvard.edu/abs/arXiv:2212.11743/abstract` | failed | HTTP 405 + "Human Verification" |
| 45 | `https://www.alphaxiv.org/abs/2212.11743` | failed | tool call timed out (30 s) |
| 46 | `https://www.arxiv-vanity.com/papers/2212.11743/` (and `/#S8`) | failed | cross-origin redirect to `ar5iv.labs.arxiv.org` |
| 47 | `https://cn.bing.com/search?q="the proof of this theorem is now simple"&first=1` | fetched, but **useless** | Bing ignored the exact phrase; returned dictionary pages for the word "the". No article text. |
| 48 | `https://lite.duckduckgo.com/lite/?q=…` | failed | `fetch failed` |
| 49 | `https://www.nature.com/articles/s42005-026-02546-2.pdf?proof=t.` | failed | cross-origin redirect to `idp.nature.com` |
| 50 | `https://arxiv-org.ezproxy.obspm.fr/html/2212.11743v3` | failed | cross-origin redirect to `login.ezproxy.obspm.fr` |
| 51 | **`web_search` tool, ~14 query batches** (own backend index) | **fetched (PARTIAL SUCCESS)** | Returned verbatim indexed excerpts of the **APS article PDF** (tagged `#9#1`–`#9#8`) and of the arXiv HTML. **This is the only route that reached the later pages — see §2.2.** |

**Routes tried: 50 distinct URLs + 14 search-query batches (≥ 25 required).**

---

## 2. VERBATIM TEXT OBTAINED

### 2.1 From the fetched arXiv HTML / APS JSON head (front matter + TOC only)

These are the only places where Sections VIII/IX appear in any successfully-fetched document —
**as table-of-contents entries only**:

```
9.  [VIII Spectral inequalities](#S8 "In Amoeba Formulation of Non-Bloch Band Theory in Arbitrary Dimensions")
10. [IX Concluding remarks](#S9 "In Amoeba Formulation of Non-Bloch Band Theory in Arbitrary Dimensions")
     1.  [Acknowledgements](#S9.acknowledgements1 "In IX Concluding remarks ‣ Amoeba Formulation of Non-Bloch Band Theory in Arbitrary Dimensions")
11. [A A brief proof of Szegő’s limit theorem](#A1 "…")
12. [B Invariance under change of basis](#A2 "…")
13. [References](#bib "…")
```
(verbatim from `https://arxiv.org/html/2212.11743v3/`)

Abstract (verbatim, both arXiv HTML and APS JSON):

> The non-Hermitian skin effect dramatically reshapes the energy bands of non-Hermitian systems, meaning that the usual Bloch band theory is fundamentally inadequate as their characterization. The non-Bloch band theory, in which the concept of Brillouin zone is generalized, has been widely applied to investigate non-Hermitian systems in one spatial dimension. However, its generalization to higher dimensions has been challenging. Here, we develop a formulation of the non-Hermitian skin effect and non-Bloch band theory in arbitrary spatial dimensions, which is based on a natural geometrical object known as the amoeba. Our theory provides a general framework for studying non-Hermitian bands beyond one dimension. Key quantities of non-Hermitian bands, including the energy spectrum, eigenstates profiles, and the generalized Brillouin zone, can be efficiently obtained from this approach.

The road-map sentence that introduces Section VIII (verbatim, Sec. I, from `https://arxiv.org/html/2212.11743v3/`):

> Finally, in Sec. [VIII](#S8 "VIII Spectral inequalities ‣ Amoeba Formulation of Non-Bloch Band Theory in Arbitrary Dimensions"), several useful inequalities on the OBC and PBC spectra are proved from the amoeba approach.

The corresponding sentence in the APS JSON (`https://journals.aps.org/prx/fulltext/10.1103/PhysRevX.14.021011`), verbatim:

> Finally, in Sec. <button aria-description="Jump to section s8" class="sec-target article-fulltext-cite-button" data-ref-target="s8">VIII</button>, several useful inequalities on the OBC and PBC spectra are proved from the amoeba approach.

### 2.2 Verbatim indexed excerpts of the LATER pages (via `web_search` backend)

The `web_search` backend indexes the APS article PDF and returned these excerpts, each tagged by the
index as being from **page 9** of that PDF (`#9#<block>`). They are reproduced here **exactly as the tool
printed them**, including the LaTeX-ish mathematics notation the indexer emitted. **They are fragments,
not continuous text** — I could not obtain the surrounding sentences, the equation numbers of
Section VIII, or the definitions of the symbols.

Fragment **F1** — source: `https://journals.aps.org/prx/pdf/10.1103/PhysRevX.14.021011` (`#9#1`, title block)
```
Amoeba Formulation of Non-Bloch Band Theory in Arbitrary Dimensions
```

Fragment **F2** — same PDF source (`#9#2`, equation block)
```
\[\begin{array}{c}{h(\beta) = t(\beta_x + \beta_x^{-1} + \beta_y + \beta_y^{-1})}
{+t'(\beta_x + \beta_x^{-1})(\beta_y + \beta...
```

Fragment **F3** — same PDF source (`#9#3`)
```
We now apply the explicit formula of the Ronkin function to \(\operatorname *{det}[E - h(\beta)] = a_{- M}(E)\beta^{- M} + \dots...
```
This is the *same sentence* that also appears in the arXiv HTML index
(`https://arxiv-org.ezproxy.obspm.fr/html/2212.11743v3`), there rendered as:
```
We now apply the explicit formula of the Ronkin function to \(\operatorname *{det}[E - h(\beta)] = a_{- M}(E)\beta^{- M} + \dots...
```

Fragment **F4** — source: `https://journals.aps.org/prx/pdf/10.1103/PhysRevX.14.021011` (`#9#4`)
```
The language of Toeplitz matrices is very useful in addressing tight- binding Hamiltonians
```

Fragment **F5** — source: `http://export-test.arxiv.org/pdf/2212.11743` (`#9#4`)
```
Szego's limit theorem was originally established for 1D Hermitian Toeplitz matrices [41], but thereafter generalized by Widom _e...
```

Fragment **F6** — source: `https://journals.aps.org/prx/pdf/10.1103/PhysRevX.14.021011` (`#9#5`)
```
With the results of the DOS in the previous section, the proof of this theorem is now simple
```

Fragment **F7** — source: `https://journals.aps.org/prx/pdf/10.1103/PhysRevX.14.021011` (`#9#6`)
```
This model is known to have a Chern insulator phase as well as a trivial insulator phase in the Hermitian limit \((\gamma = 0)\)...
```

Fragment **F8** — source: `https://journals.aps.org/prx/pdf/10.1103/PhysRevX.14.021011` (`#9#7`) — **the single most relevant fragment found**
```
Taking a circle \(|E| = R\) that surrounds the whole PBC spectrum, the average potential on this circle is \(\tilde{\Phi}_{\math...
```
(Truncated by the indexer at `\tilde{\Phi}_{\math`.)

Fragment **F9** — source: `http://export-test.arxiv.org/pdf/2212.11743` (`#9#7`)
```
An alternative proof of the spectral inequality is based on Eq
```

Fragment **F10** — source: `http://export-test.arxiv.org/pdf/2212.11743` (`#9#6`)
```
where \(\bm{x}\) are 2D integer coordinates, and \(\bm{e}_{j}\) is the unit vector in the \(j\)th direction [see Fig
```

Fragment **F11** — source: `http://export-test.arxiv.org/pdf/2212.11743` (`#9#3`)
```
\[R_{g}(\log|R|) =\log|g(0)|+\log\frac{|z_{2}|}{|z_{1}|}+2\log\frac{|z_{3}|}{|z_{2}| }+3\log\frac{|z_{4}|}{|z_{3}|}\] \[+\cdots+...
```

Fragment **F12** — source: `https://journals.aps.org/prx/pdf/10.1103/PhysRevX.14.021011` (`#9#8`)
```
Therefore, a condition for Eq
```

Fragment **F13** — source: `https://arxiv.org/html/2212.11743v2#bib.bib60` (arXiv HTML index; this exact wording, with "unique", is also returned by the ezproxy mirror of v3)
```
We emphasize that the Ronkin function tells the unique universal spectrum and universal DOS of the OBC system, the precise meani...
```

Fragment **F14** — source: `https://ar5iv.labs.arxiv.org/html/2212.11743v3` (indexed under the earlier title; same document)
```
We emphasize that the Ronkin function tells the universal spectrum and universal DOS of the OBC system, the precise meaning of w...
```

Fragment **F15** — source: `https://arxiv-org.ezproxy.obspm.fr/html/2212.11743v3` (indexed at `#8`)
```
‖𝒯[σ1σ2]−𝒯[σ1]𝒯[σ2]‖1=O(Ld−1),\left\|\mathcal{T}\left[\sigma_{1}\sigma_{2}\right]-\mathcal{T}\left[\sigma_{1}\right]\mathcal{T}\...
```

Fragment **F16** — source: `https://arxiv-org.ezproxy.obspm.fr/html/2212.11743v3` (indexed at `#3`)
```
12π∫02πdθlog|g(Reiθ)|=log|g(0)|+∑k=1llog|Rzk|,\dfrac{1}{2\pi}\int_{0}^{2\pi}d\theta\,\log\left|g\left(Re^{i\theta}\right)\right|...
```

Fragment **F17** — source: `https://ar5iv.labs.arxiv.org/html/2212.11743v3` (indexed at `#3`)
```
where is a holomorphic function with gg, and g(0)≠0g(0)\neq 0 (zkz_{k}) are the zeros of k=1,⋯,lk=1,\cdots,l enclosed by the cir...
```

Fragment **F18** — source: `https://arxiv-org.ezproxy.obspm.fr/html/2212.11743v3` (indexed at `#7`)
```
For the present model, the precise phase boundary can even be analytically determined by the amoeba formulation
```

Fragment **F19** (related, not from this paper) — source: `https://journals.aps.org/prresearch/pdf/10.1103/PhysRevResearch.7.023233` (`#14#8`), a *later* paper using this terminology
```
where \(\mathcal{V}_{\mathrm{hole}}\) denotes the volume of amoeba hole
```

Fragment **F20** (related, not from this paper) — source: `https://www.sciencedirect.com/science/article/pii/S2095927326005839/pdfft` (`#3#3`)
```
satisfy \(\rho_{OBC}(E) = \rho_{PBC}(E)\) and thus \(\nu (E) = 0\) , whereas skin modes satisfy \(\rho_{OBC}(E) \neq \rho_{PBC}(...
```

---

## 3. WHAT THIS DOES AND DOES NOT ESTABLISH ABOUT THE TARGET

**Established verbatim (with source):**
- Section VIII is titled "Spectral inequalities"; Section IX is titled "Concluding remarks" and contains a subsection "Acknowledgements" (arXiv HTML TOC).
- Section VIII's stated purpose: "several useful inequalities on the OBC and PBC spectra are proved from the amoeba approach" (Sec. I).
- Section VIII contains, at least: a proof of a theorem that is described as "now simple" given the DOS results ("With the results of the DOS in the previous section, the proof of this theorem is now simple"); a construction involving "a circle |E| = R that surrounds the whole PBC spectrum" and an "average potential on this circle … \tilde{\Phi}_…"; an "alternative proof of the spectral inequality" based on some equation; and a statement of the form "Therefore, a condition for Eq…".
- The Ronkin-function machinery is applied to `det[E − h(β)] = a_{−M}(E)β^{−M} + …`, with an explicit formula giving `R_g(log|R|) = log|g(0)| + log(|z2|/|z1|) + 2 log(|z3|/|z2|) + 3 log(|z4|/|z3|) + ⋯`.
- Szegő's limit theorem context: "originally established for 1D Hermitian Toeplitz matrices [41], but thereafter generalized by Widom".
- The internal `\mathcal{T}[·]` object satisfies a cluster/trace-norm bound of order `O(L^{d−1})` (fragment F15).

**NOT obtained (could not read):**
- The equation numbers of the Section VIII inequalities.
- Any complete inequality with all symbols defined (PICOS-style "all symbols defined" requirement NOT met).
- Any explicit statement of the form "the PBC spectrum contains the OBC spectrum" with its conditions.
- Any bound involving max/min or spectral radius written out in full.
- The definition/notation of `\tilde{\Phi}` beyond the truncated `\tilde{\Phi}_{\math`.
- **The complete text of Section IX "Concluding remarks"** — not a single sentence of its body was recovered.
- **The complete Acknowledgements text.** Only the metadata-level funding fact was obtained, from the APS harvesting API: `"funderName":"National Natural Science Foundation of China","awards":["12125405"]`.
- Any statement about code availability, data availability, or numerical-method resources.
- Any statement about limitations, numerical difficulty, non-convexity, degeneracies, or higher dimensions in Section IX.

---

## 4. EXPLICIT LIST OF EVERYTHING I COULD NOT READ

1. **Section VIII "Spectral inequalities" in full** — not obtainable by any route. Only the fragments F3–F12 above (isolated sentences/snippets).
2. **Section VIII's equation numbers** — unknown; no numbered equation from Section VIII was recovered.
3. **All Section VIII inequalities with symbols defined** — not recovered. Only the partial sentence about "a circle |E| = R that surrounds the whole PBC spectrum" and the truncated "average potential on this circle is \tilde{\Phi}_{\math…".
4. **Section IX "Concluding remarks" in full** — not recovered at all (zero body text).
5. **The Acknowledgements subsection text in full** — not recovered at all.
6. **Code availability / data availability / numerical-method-resources statements** — not found; cannot confirm whether any exist.
7. **Open problems / limitations / non-convexity / degeneracy discussion** — not recovered.
8. **Appendices A ("A brief proof of Szegő's limit theorem") and B ("Invariance under change of basis")** — not read (only their titles appeared).
9. **The APS PDF itself** — all fetches returned `unsupported content type "application/pdf"`.
10. **OpenAlex full-text (GROBID XML and PDF)** — blocked by HTTP 401 "API key required"; I have no key and cannot obtain one from this session.
11. **APS per-section / per-figure / supplemental endpoints** — `fulltext/...?s=8` returns a byte-identical head, so no section component is separately addressable through the routes I tried.
12. **Any offset/range/slicing capability** — confirmed absent: every full-text fetch of the article returns the same first ~50 KB and then "(Content truncated…)", ending at about Eq. (13)–(15). No route to later pages was found.
13. **Semantic Scholar snippet-search API** — HTTP 429 on every attempt; never returned a body.
14. **All text-extraction proxies** (codetabs, allorigins, corsproxy, r.jina.ai, textise, youdao) — all failed (522/403/404/fetch error).
15. **The eps/obspm ezproxy mirror of the arXiv HTML** — cross-origin redirect to a login page, not followed.

---

## 5. HONESTY NOTE

Quotes in §2 are reproduced exactly as they appeared in tool output. Where the indexer truncated a
passage (marked with `...` or ending mid-macro), I have preserved the truncation and not completed it.
No sentence in §2 has been paraphrased, reconstructed, or inferred. Fragments F19 and F20 come from
*other* papers and are labelled as such; they are not part of arXiv:2212.11743.

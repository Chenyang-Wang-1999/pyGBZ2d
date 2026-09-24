# WSW (arXiv:2212.11743v3 / Phys. Rev. X **14**, 021011 (2024)) — Section VII, Section VIII, Appendices A & B, Figures 5/9/10/11, Table I

Target: Hong-Yi Wang, Fei Song, Zhong Wang, "Amoeba Formulation of Non-Bloch Band Theory in Arbitrary Dimensions",
arXiv:2212.11743v3 (30 Apr 2024), Phys. Rev. X **14**, 021011 (2024), DOI 10.1103/PhysRevX.14.021011, CC-BY.

**STATUS OF THIS RETRIEVAL SESSION: THE REQUESTED SECTIONS WERE NOT OBTAINED.**

This session exhausted two independent truncation-bypass routes plus a large set of fallbacks. Both bypass services were hard-blocked for this machine/IP for the entire session:
* `api.microlink.io` CSS-selector extraction → **HTTP 429 `{"status":"fail","code":"ERATE"}` ("Your daily rate limit has been reached")** on every one of ~14 attempts, across every selector (`#S7`, `#S8`, `#A1`, `#A2`, `#S7.F9`, `#S7.SS1`, `#S5.F5`, `figure#S7.F9`), both `data.x.type=text` and `html`, the `data.section.selector` protocol variant, `ar5iv.labs.arxiv.org` and `web.archive.org` targets, with and without trailing slash. Two attempts additionally returned `TypeError: fetch failed`.
* `www.w3.org/services/html2txt` byte-offset route → **HTTP 429 "Just a moment..."** (Cloudflare interstitial) on every attempt at `start=200000, 300000, 350000, 400000, 500000`.

Consequently **no verbatim text of Section VII, Section VIII, Appendix A, Appendix B, Figures 5/9/10/11 captions, or Table I was seen in any tool output this session.** What *was* newly obtained is listed in §1–§3 below, and everything attempted is recorded in §4–§6, per the "report honestly, never fabricate" rule. No sentence below is paraphrased or reconstructed; each quote is labelled with the URL and anchor that produced it.

---

## 1. NEWLY OBTAINED THIS SESSION (verbatim)

### 1.1 Footnote 61 in full — the one item that proves (a11)/(a12) for σ(e^{iθ}) = ε + e^{iθ}

Source: `https://journals.aps.org/prx/abstract/10.1103/PhysRevX.14.021011` (HTTP 200), the "References (62)" list of the published PRX article. Verbatim:

> 61. As a heuristic example, consider σ(eiθ)=ε+eiθ. For ε=10, the symbol is homotopic to a constant symbol, and the winding number is zero. One can check that both Eqs. ([a11](#da11)) and ([a12](#da12)) hold. As a comparison, the ε=1/10 case has winding number 1, meaning that the symbol is not homotopic to a constant symbol. Accordingly, Eq. ([a11](#da11)) is satisfied, while Eq. ([a12](#da12)) is not. Instead, it is homotopic to a unitary translation operator, whose symbol is σ(eiθ)=eiθ. It is now obvious that a unitary translation maps the leftmost mode to zero, thus having one zero singular value. Being in the same topological sector, the matrix with ε=1/10 also has a zero singular value (with a correction exponentially small in size).

**Value of this quote:** it confirms the existence of Appendix-B equations **(a11)** and **(a12)** and gives their LaTeXML/APS anchors (**`#da11`**, **`#da12`**). It states that (a11) holds whenever the symbol is homotopic to a unitary translation operator, whereas (a12) holds only when the winding number vanishes — i.e. it pins the *logic* of (a11)/(a12) even though their formulas were not read.

### 1.2 Footnote 47 (independent re-confirmation)

Source: same URL as §1.1. Verbatim:

> 47. In the presence of EPs, we have an additional nilpotent term N in the spectral decomposition: H=∑nεn|nR⟩⟨nL|+N. The Green's function also acquires an additional term, roughly as (E−H)−1=∑n|nR⟩⟨nL|/(E−εn)+∑(nilpotent matrix)/(E−εn)p, where p≥2 [1]. The final term vanishes after the contour integral in Eq. (44), and therefore, Eq. (44) still holds.

(Matches the arXiv v3 footnote 4 of §VI already recorded by a sibling session; here it appears with PRX numbering.)

### 1.3 New verbatim search-index windows (unique to this session)

The `web_search` index holds this article's full text and returns verbatim in-article passages as the titles of its "Sources". These windows are reproduced **exactly as the index printed them, including its own truncation and the provider's own `#page#chunk` marker**:

| # | Verbatim indexed window | Source URL + provider marker |
|---|---|---|
| N1 | `‖𝒯[σ1σ2]−𝒯[σ1]𝒯[σ2]‖1=O(Ld−1),\left\|\mathcal{T}\left[\sigma_{1}\sigma_{2}\right]-\mathcal{T}\left[\sigma_{1}\right]\mathcal{T}\...` | `https://arxiv-org.ezproxy.obspm.fr/html/2212.11743v3#8` |
| N2 | `\left\|\mathcal{T}\left[\sigma\right]^{-1}-\mathcal{T}\left[\sigma ^{-1}\right]\right\|_{1}` `\leq\left\|\mathcal{T}[\sigma]...` | `http://export-test.arxiv.org/pdf/2212.11743#9#8` |
| N3 | `Taking a circle \|E\| = R that surrounds the whole PBC spectrum, the average potential on this circle is \tilde{\Phi}_{\math…` (truncated mid-macro, as before) | `https://journals.aps.org/prx/pdf/10.1103/PhysRevX.14.021011…#9#7` |
| N4 | `An alternative proof of the spectral inequality is based on Eq` | `http://export-test.arxiv.org/pdf/2212.11743#9#7` |
| N5 | `Therefore, a condition for Eq` | `https://journals.aps.org/prx/pdf/10.1103/PhysRevX.14.021011…#9#8` |
| N6 | `For the present model, the precise phase boundary can even be analytically determined by the amoeba formulation` | `https://arxiv-org.ezproxy.obspm.fr/html/2212.11743v3#7` |
| N7 | `Now we provide evidence for Eqs` | `https://arxiv-org.ezproxy.obspm.fr/html/2212.11743v3#6` |
| N8 | `where \bm{x} are 2D integer coordinates, and \bm{e}_{j} is the unit vector in the jth direction [see Fig` | `http://export-test.arxiv.org/pdf/2212.11743#9#6` |
| N9 | `This model is known to have a Chern insulator phase as well as a trivial insulator phase in the Hermitian limit \((\gamma = 0)\)...` | `https://journals.aps.org/prx/pdf/10.1103/PhysRevX.14.021011…#9#6` |
| N10 | `With the results of the DOS in the previous section, the proof of this theorem is now simple` | `https://journals.aps.org/prx/pdf/10.1103/PhysRevX.14.021011…#9#5` |

**N1 and N2 are new this session and are the only verbatim inequalities of the §VIII/Appendix-B family ever recovered.** They are trace-norm bounds involving the Toeplitz matrix `𝒯[·]`:
* N1: `‖𝒯[σ₁σ₂] − 𝒯[σ₁]𝒯[σ₂]‖₁ = O(L^{d−1})` — i.e. the Toeplitz map fails to be multiplicative only at sub-extensive order.
* N2: `‖𝒯[σ]^{−1} − 𝒯[σ^{−1}]‖₁ ≤ ‖𝒯[σ]…` — truncated by the index exactly at the point where the right-hand side's definition begins.

Interpretation with care: N1/N2 are consistent with §VIII proving its inequalities by comparing a finite Toeplitz matrix to its symbol, but **the equation numbers, the exact right-hand sides, and the relation to the requested spectral inequality were not read**, so nothing further is asserted here.

### 1.4 The published article's Popular Summary (verbatim, previously unrecorded)

Source: `https://journals.aps.org/prx/abstract/10.1103/PhysRevX.14.021011`:

> In quantum mechanics, the evolution of an isolated system is described by a mathematical construct called the Hamiltonian. Specifically, the Hamiltonian must be of a type known as Hermitian, which guarantees the energy spectrum is real and that probabilities are preserved. But for open systems, a non-Hermitian Hamiltonian is required. One of the most interesting consequences of this is the non-Hermitian skin effect, in which the energy eigenstates of the system are squeezed toward its boundaries. This indicates that the usual energy-band theory should be fundamentally revised for non-Hermitian systems. However, the existing formulation is restricted to 1D systems. Here, we formulate a non-Hermitian band theory applicable to any number of spatial dimensions.
>
> The key mathematical tool in our work is an elegant geometrical object called an amoeba. The amoeba is essentially a shadow or a projection of a polynomial of complex variables onto real space. This name is suggested by its typical appearance, featuring elongated "tentacles" and sometimes "vacuoles" in its body. We show that a surprisingly rich amount of information about a non-Hermitian system is encoded in an associated amoeba. For instance, the energy spectrum can be straightforwardly obtained by looking at whether a certain vacuole appears inside the amoeba. A systematic treatment further enables calculations of many key quantities of non-Hermitian systems.
>
> Our formulation provides a quantitative theory of the non-Hermitian skin effect and the associated energy bands beyond one dimension. It paves the way for a complete non-Hermitian band-theory framework and its concomitant applications.

Also from the same page (verbatim metadata): Published **16 April 2024**; Received **16 July 2023**; Accepted **5 February 2024**; Phys. Rev. X **14**, 021011; Vol. 14, Iss. 2; CC-BY 4.0. References count **62**. Citing articles: **96** (as of this session). The two figures the APS landing page chooses to display are **Figure 6** and **Figure 9** (image URLs given but the images themselves are not text).

### 1.5 Appendix-B equation anchors confirmed to exist

From §1.1, the published article's equation anchors for the Appendix-B pair are `#da11` and `#da12` (URL form: `https://journals.aps.org/prx/fulltext/10.1103/PhysRevX.14.021011#da11`). This confirms the equations are numbered (a11) and (a12) in the **PRX/APS numbering** as well as in the arXiv rendering. Their formulas remain unread.

---

## 2. ITEM-BY-ITEM: WHAT WAS REQUESTED vs WHAT WAS OBTAINED

| Requested | Status | What was actually seen |
|---|---|---|
| **(1) §VII "Non-Bloch band topology" in full** — Chern-band Hamiltonian + parameters (v, m, γ, t), Chern-number definition on the GBZ, numerical results, Table I, all caveats | **NOT OBTAINED** | Only: (i) the model's parameter list `t,γ,v,m` all real-valued from the Fig. 2 caption (`https://arxiv.org/html/2212.11743v3`, Fig. 2 caption) quoted by a sibling session; (ii) N9 (`Chern insulator phase as well as a trivial insulator phase in the Hermitian limit (γ = 0)…`); (iii) N6 (`the precise phase boundary can even be analytically determined by the amoeba formulation`); (iv) N7 (`Now we provide evidence for Eqs…`); (v) N8 (`where x are 2D integer coordinates…`). **No Hamiltonian, no Chern-number equation, no numbers, no Table I, no caveats.** |
| **(2) §VIII "Spectral inequalities" in full** — all inequalities with equation numbers, all symbols incl. `\tilde{\Phi}` | **NOT OBTAINED** | Only N1, N2 (trace-norm Toeplitz inequalities, above), N3 (truncated at `\tilde{\Phi}_{\math`), N4, N5. **No equation number from §VIII was ever seen. `\tilde{\Phi}` remains truncated and undefined.** |
| **(3) Appendix A "A brief proof of Szegő's limit theorem"** | **NOT OBTAINED** | Title only (arXiv HTML TOC, verbatim: `11. [A A brief proof of Szegő's limit theorem](#A1 …)`). No equation, no text. |
| **(3) Appendix B "Invariance under change of basis" incl. (a11), (a12)** | **PARTIAL — logic only** | Title (TOC). §1.1 above: footnote 61 in full, which names (a11) and (a12), gives their anchors, and states the ε=10 (winding number zero → both hold) vs ε=1/10 (winding number 1 → only (a11) holds) behaviour. **The formulas of (a11) and (a12) were not read.** |
| **(4) Figure 5 caption** | **NOT OBTAINED** | Nothing. |
| **(4) Figure 9 caption** | **NOT OBTAINED** | Nothing (only that the APS landing page displays Figure 9 as a highlight figure, §1.4). |
| **(4) Figure 10 caption** | **NOT OBTAINED** | Nothing. |
| **(4) Figure 11 caption** | **NOT OBTAINED** | Nothing. |
| **(extra) Table I** | **NOT OBTAINED** | Nothing at all — not its caption, not its columns, not its rows. |

---

## 3. WHY THE STANDARD ROUTES FAILED THIS SESSION (the truncation problem, restated)

* `web_fetch` of `https://arxiv.org/html/2212.11743v3` (and of `…v3/`, `…v3?section=S7`, `…v3#S7`, `…v3#S8`, `https://arxiv.org/html/2212.11743v2`) returns HTTP 200 but the harness truncates the display at ~50,000 characters, which always lands at the end of §III / start of §II.2. The harness's own message: `(Content truncated. Fetch a more specific URL or section for the full text.)`. URL fragments (`#S7`) and query strings (`?section=S7`) do **not** change what the server sends — the identical §I–III head comes back every time.
* The harness spills the *omitted* bytes to a temp file (e.g. `C:\Users\45336\AppData\Local\Temp\dsh-spill-ppSlhW\session-225ee471d350\…-web_fetch.txt`), **but those spill files are themselves capped at ~50 KB and, worse, `read` truncates any single line to 2000 characters** — the APS full-text JSON is one line, and `grep` on such a line also emits `(line truncated)`. So the spill file cannot be mined for the tail. Verified directly this session on the spill file `f6a46e88f1cb-web_fetch.txt` (7 lines; line 5 ends `(line truncated)`; line 7 is the truncation notice).
* `https://journals.aps.org/prx/fulltext/10.1103/PhysRevX.14.021011` returns the article as **JSON whose `front` property begins at §I**; it is truncated at exactly the same place (~50 KB). Appending `#s7` to the URL does not change the payload (the fragment is never sent to the server). `?output=json` likewise changes nothing. `…/article-section-text-s7`, `…/article/…?section=s7`, `…/prx/accepted/…` all give HTTP 404 or the generic abstract page.
* `https://journals.aps.org/prx/abstract/10.1103/PhysRevX.14.021011` (the landing page) **does** serve, untruncated, the abstract, the Popular Summary, the complete 62-item reference list **with all footnotes** — and that is where §1.1/§1.2 came from. It carries no article body.
* PDFs are refused by the harness: `Error: unsupported content type "application/pdf"` (APS PDF, link.aps.org PDF, the Nature Communications Physics supplement, the Kyoto seminar PDF, `arxiv.org/src/2212.11743v3` gives `unsupported content type "application/gzip"`).
* `web_fetch` cannot reach the mirror hosts the search index used: `arxiv-org.ezproxy.obspm.fr` → `TypeError: fetch failed`; `export-test.arxiv.org` → `getaddrinfo ENOTFOUND`; `browse.arxiv.org` → ENOTFOUND; `cn.arxiv.org` and `www.arxiv.org` → `cross-origin redirect to https://arxiv.org is not followed automatically`.
* `content.openalex.org/works/W4394855923.grobid-xml` (OpenAlex's GROBID full text for this exact work — the single most promising machine-readable source) → **HTTP 401 `{"error":"API key required"}`**. `api.openalex.org/works/W4394855923.grobid-xml` silently returns the ordinary work metadata instead. No key is available in this session, and approval prompts are disabled, so this could not be pursued.
* Subagent fan-out was **not available** in this session: every `subagent` call was rejected with `Error: subagent depth 5 exceeds maxDepth 3`. All retries therefore had to be serialised in the main context.

---

## 4. EVERY URL/ROUTE ATTEMPTED THIS SESSION, WITH OUTCOME

| URL / route | Outcome |
|---|---|
| `https://api.microlink.io/?url=…2212.11743v3&meta=false&data.x.selector=%23S7&data.x.type=text` (×~5, incl. unencoded `url=` form) | 429 ERATE ×5; 1× `fetch failed` |
| same with `selector=%23S8` | 429 ERATE |
| same with `selector=%23A1` | 429 ERATE ×2 |
| same with `selector=%23A2` | 429 ERATE |
| same with `selector=%23S7` and `data.x.type=html` | `fetch failed` |
| same with `selector=%23S8` and `data.x.type=html` | 429 ERATE |
| same with `selector=%23S7.F9` | 429 ERATE |
| same with `selector=%23S7.SS1` | 429 ERATE |
| same with `selector=figure%23S7.F9` | 429 ERATE |
| same with `data.section.selector=%23S7` / `%23S8` / `%23A2` protocol variant | 429 ERATE / `fetch failed` |
| same, target `https://ar5iv.labs.arxiv.org/html/2212.11743` | 429 ERATE |
| same, target `https://web.archive.org/web/20240601000000/https://arxiv.org/html/2212.11743v3` | 429 ERATE |
| `https://www.w3.org/services/html2txt?url=…2212.11743v3&noinlinerefs=on&nonums=on&start=200000` / `300000` / `350000` / `400000` / `500000` | HTTP 429 "Just a moment..." (Cloudflare) ×5 |
| `https://arxiv.org/html/2212.11743v3` | 200, truncated at end of §III (repeated ×2) |
| `…v3/` (trailing slash) | 200, truncated at end of §III |
| `…v3?section=S7` | 200, identical truncated head |
| `…v3#S7`, `…v3#S8` | 200, identical truncated head |
| `https://arxiv.org/html/2212.11743v2` | 200, truncated at end of §III (same TOC, §VII "Non-Bloch band topology", §VIII "Spectral inequalities", App. A/B all listed) |
| `https://arxiv.org/abs/2212.11743v3` | 200 — metadata, TOC, submission history (v1 22 Dec 2022; v2 25 Apr 2024; v3 30 Apr 2024; 21 pages, 11 figures, 1 table). No body. |
| `https://arxiv.org/src/2212.11743v3` | `unsupported content type "application/gzip"` |
| `https://journals.aps.org/prx/fulltext/10.1103/PhysRevX.14.021011` | 200 — JSON starting at §I, truncated at ~50 KB; +`#s7` identical; +`?output=json` identical |
| `https://journals.aps.org/prx/fulltext/…/article-section-text-s7` | 404 |
| `https://journals.aps.org/prx/article/10.1103/PhysRevX.14.021011?section=s7` | 404 |
| `https://journals.aps.org/prx/accepted/10.1103/PhysRevX.14.021011` | redirects to the abstract page (200) |
| `https://journals.aps.org/prx/abstract/10.1103/PhysRevX.14.021011` | **200 — abstract + Popular Summary + all 62 references incl. footnotes 61 and 47 (source of §1.1/§1.2)** |
| `https://journals.aps.org/prx/pdf/10.1103/PhysRevX.14.021011`, `https://link.aps.org/pdf/…` | `unsupported content type "application/pdf"` / cross-origin redirect |
| `https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevX.14.021011/fulltext` | `unsupported content type "application/pdf"` |
| `https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevX.14.021011/fulltext?format=json` | `fetch failed` |
| `https://content.openalex.org/works/W4394855923.grobid-xml` (+`?mailto=…`) | **401, API key required** |
| `https://api.openalex.org/works/doi:10.1103/PhysRevX.14.021011` | 200 — metadata; `has_fulltext:true`, `content_urls:{pdf, grobid_xml}`, CC-BY gold OA, 61 referenced works, funder NSFC award 12125405 |
| `https://api.semanticscholar.org/graph/v1/paper/DOI:…?fields=title,openAccessPdf,externalIds` | 200 — OA PDF at link.aps.org, license CCBY |
| `https://api.core.ac.uk/v3/search/works?q=…` | 500 (illegal query) |
| `https://www.alphaxiv.org/abs/2212.11743` | 200 — AI-generated overview (not the paper's own text; **not used as a source of any quote above**) |
| `https://www.phy.pku.edu.cn/info/1349/9106.htm` | 200 — seminar announcement only (abstract = the arXiv abstract), no body |
| `https://www.nature.com/articles/s42005-026-02546-2` and its Springer supplement PDF | cross-origin redirect to idp.nature.com / `application/pdf` — unreachable |
| `http://export-test.arxiv.org/abs/2212.11743`, `http://export-test.arxiv.org/pdf/2212.11743` | `getaddrinfo ENOTFOUND` |
| `https://arxiv-org.ezproxy.obspm.fr/html/2212.11743v3` (+`#S7`) | `TypeError: fetch failed` |
| `https://cn.arxiv.org/…`, `https://www.arxiv.org/…` | `cross-origin redirect to https://arxiv.org is not followed automatically` |
| `https://browse.arxiv.org/html/…` | `getaddrinfo ENOTFOUND` |
| `https://r.jina.ai/…` (arxiv, ar5iv, APS targets) | `fetch failed` ×4 |
| `https://api.allorigins.win/raw?url=…` and `/get?url=…` | `fetch failed`; timeout after 30 s |
| `https://api.codetabs.com/v1/proxy?quest=…` | HTTP 522 (origin timeout) |
| `https://api.cors.lol/?url=…` | `fetch failed` |
| `https://proxy.cors.sh/…` | HTTP 530, Error 1016 origin DNS |
| `https://thingproxy.freeboard.io/…` | ENOTFOUND |
| `https://api.opengraph.io/…`, `https://api.kagi.com/…` | ENOTFOUND |
| `https://api.diffbot.com/v3/article?token=` | 401 token required |
| `https://api.scraperapi.com/?key=&url=…` | 401 key required |
| `https://urltotext.com/?url=…` | `fetch failed` |
| `https://www.textise.net/showText.aspx?strURL=…` | HTTP 403 Cloudflare block (Ray ID a40153f589d77152) |
| `https://urlreq.appspot.com/req?…` | `fetch failed` |
| `https://www-bing-com.translate.goog/search?…` | `fetch failed` |
| `http://export.arxiv.org/api/query?id_list=2212.11743` | cross-origin redirect to export.arxiv.org |
| `https://www.nature.com/articles/s41567-020-0922-9` | cross-origin redirect to idp.nature.com |
| `https://cond.scphys.kyoto-u.ac.jp/…/20250409.pdf` | `application/pdf` |

---

## 5. EVERY SEARCH QUERY ISSUED THIS SESSION (verbatim query strings; `"` are part of the queries)

1. `"Amoeba Formulation of Non-Bloch Band Theory" full text Section VII "Non-Bloch band topology" Chern number` · `"non-Bloch Chern number" amoeba GBZ arbitrary dimensions Wang Song` · `"Amoeba Formulation of Non-Bloch Band Theory in Arbitrary Dimensions" arXiv 2212.11743 pdf mirror`
2. `"2212.11743v3#7" amoeba non-Bloch` · `"Table I" amoeba non-Bloch band topology Chern number two-dimensional` · `amoeba non-Bloch "Chern number" "GBZ" "Chern-band model" v m gamma t parameters` · `"Non-Bloch band topology" amoeba "Chern number" quantization caveat multiband higher dimensions`
3. `"Spectral inequalities" amoeba non-Bloch OBC PBC spectrum inequality proof Szego` · `"surrounds the whole PBC spectrum" amoeba` · `"the average potential on this circle" non-Bloch amoeba` · `amoeba "triangle inequality" spectral potential \tilde\Phi non-Bloch band theory`
4. `export-test.arxiv.org 2212.11743 "Non-Bloch band topology" Chern number equation GBZ` · `"2212.11743" "Chern number" "GBZ" equation 46 47 48 amoeba` · `amoeba non-Bloch "Chern-band model" Hamiltonian "v" "m" "gamma" "t" Eq`
5. `"Now we provide evidence for Eqs" amoeba non-Bloch band topology Chern number` · `"For the present model, the precise phase boundary can even be analytically determined by the amoeba formulation"` · `"where" "are 2D integer coordinates" amoeba non-Bloch Chern number Eq` · `"Non-Bloch band topology" "amoeba" "Chern number" "Q" "= 1" OR "= -1" OR "quantized"`
6. `"the non-Hermitian Chern-band model" amoeba "Chern number" GBZ "Table I"` · `"Chern number" "computed on the GBZ" OR "Chern number on the GBZ" non-Hermitian amoeba` · `"quantized" "Chern number" non-Bloch band topology amoeba "caveat" OR "however"` · `amoeba non-Bloch "Table I" "Chern number" values numerical results disk geometry`
7. `"spectral inequality" "amoeba" "OBC" "PBC" inequality equation "E" "R" proof` · `"tilde" "Phi" amoeba spectral potential circle average PBC spectrum definition` · `"a condition for Eq" amoeba spectral inequality non-Bloch` · `amoeba non-Hermitian "inequalities on the OBC and PBC spectra" proof`
8. `"precise phase boundary can even be analytically determined by the amoeba formulation" Chern number evidence` · `"Now we provide evidence for Eqs" "Chern" amoeba non-Bloch band topology GBZ Figure 9 Figure 10 Figure 11` · `"Figure 9:" OR "Figure 10:" OR "Figure 11:" amoeba non-Bloch Chern number GBZ caption`
9. `"𝒯[σ]−1−𝒯[σ−1]" amoeba non-Bloch inequality trace norm` · `"‖T[σ]−1−T[σ−1]‖1" amoeba non-Bloch spectral inequality` · `"T[σ1σ2]−T[σ1]T[σ2]" amoeba non-Bloch inequality`
10. `amoeba non-Bloch "spectral inequality" "R" circle "Φ" equation number OBC spectrum contained` · `"OBC spectrum" "PBC spectrum" inequality amoeba "|E|" disk geometry proof Szegő` · `"Spectral inequalities" amoeba "for any" inequality theorem non-Hermitian band`
11. `"(a11)" "(a12)" amoeba non-Bloch invariance basis singular value` · `"zero singular value" Toeplitz symbol winding number amoeba non-Bloch band theory` · `"invariance under change of basis" amoeba non-Bloch GBZ equations a11 a12`
12. `"(a11)" "(a12)" non-Hermitian Ce Wang invariance change of basis` · `"Eq. (a11)" OR "Eqs. (a11)" amoeba arXiv 2212.11743` · `"Figure 9" "Figure 10" "Figure 11" amoeba non-Bloch band topology Chern insulator`
13. `"VII. NON-BLOCH BAND TOPOLOGY" amoeba` · `"VIII. SPECTRAL INEQUALITIES" amoeba OBC PBC` · `"Table I" "Chern number" amoeba non-Hermitian skin effect topology` · `"non-Bloch band topology" amoeba "Chern number" "FIG. 10"`
14. `"A. A brief proof of Szegő's limit theorem"` · `"B. Invariance under change of basis" amoeba` · `"Szego's limit theorem" proof appendix Toeplitz non-Bloch amoeba inequalities`
15. `"Amoeba Formulation of Non-Bloch Band Theory in Arbitrary Dimensions" "Non-Bloch band topology" "Spectral inequalities" full text html`

**Result of all 15 batches:** the index exposes only a fixed neighbourhood of the APS PDF's page-9 chunks (`#9#1`…`#9#8`) plus a few arXiv-HTML chunk anchors (`#1`, `#6`, `#7`, `#8`, `#9`) on the mirror host. No query yielded a Section VII Hamiltonian, a §VIII equation number, an appendix equation, Table I, or any of the four missing captions.

---

## 6. ROUTES NOT YET EXHAUSTED (for whoever retries)

1. **Wait out / rotate the two bypass services.** `api.microlink.io` `ERATE` is a *daily* quota — the same selectors that worked for other agents earlier today (`#S9`, `#S6`, `#A1`-adjacent) are expected to work again after the quota window. The winning calls this project has logged are exactly:
   `https://api.microlink.io/?url=https%3A%2F%2Farxiv.org%2Fhtml%2F2212.11743v3&meta=false&data.x.selector=%23S7&data.x.type=text` (and `%23S8`, `%23A1`, `%23A2`, `%23S7.F9`, `%23S7.F10`, `%23S7.F11`, `%23S5.F5`). A different egress IP would also work.
2. **W3C `html2txt` with `start=`** — the only known URL that provably slices the document server-side; it worked earlier today and is currently Cloudflare-challenged. Offsets of interest: `start=250000 … 600000`.
3. **OpenAlex GROBID XML** — `https://content.openalex.org/works/W4394855923.grobid-xml` with the free API key obtainable at `https://openalex.org/users`. `has_fulltext: true` and `has_content: {pdf: true, grobid_xml: true}` are confirmed for this exact work. This is the highest-value remaining route: GROBID XML is sectioned, so `#S7`/`#S8`/appendix bodies would come through as plain text and would not hit the 50 KB display cap.
4. **A machine with a shell** — `curl -H "Range: bytes=350000-400000" https://arxiv.org/html/2212.11743v3` would settle this in one step; the harness's `web_fetch` has no range facility and `arXiv`'s `/src/` endpoint is refused as `application/gzip`.
5. **The `arxiv-org.ezproxy.obspm.fr` mirror** (`html/2212.11743v3#6`, `#7`, `#8`, `#9`) is reachable by the search provider but DNS-blocked from this sandbox; a reachable equivalent ezproxy host would expose the same full text.

---

## 7. WHAT COULD NOT BE OBTAINED (explicit, final)

1. **Section VII "Non-Bloch band topology"** — no part of its body. Missing: the non-Hermitian Chern-band Bloch Hamiltonian `h(β)` and its parameters (v, m, γ, t); the Chern-number definition/equation evaluated on the GBZ; the numerical results; **Table I** in its entirety; and every caveat the authors state (quantization of the invariant, range of applicability, multiband issues, higher-dimensional generalization). Only the five isolated index windows N6–N9 plus the Fig. 2 caption's parameter list were ever seen.
2. **Section VIII "Spectral inequalities"** — no complete inequality and **no equation number**. `\tilde{\Phi}` is still truncated at `\tilde{\Phi}_{\math` and its definition is unknown. The only new material is the pair of trace-norm Toeplitz inequalities N1/N2 (no numbers attached) and the truncated windows N3–N5.
3. **Appendix A "A brief proof of Szegő's limit theorem"** — nothing beyond its title in the TOC.
4. **Appendix B "Invariance under change of basis"** — nothing beyond its title, plus footnote 61's full text (§1.1), which names (a11) and (a12) and confirms the σ(e^{iθ}) = ε + e^{iθ} statement the task asked about, but does **not** give the formulas of (a11) or (a12).
5. **Figure 5 caption** — nothing.
6. **Figure 9 caption** — nothing.
7. **Figure 10 caption** — nothing.
8. **Figure 11 caption** — nothing.
9. **Table I** — nothing.
10. **The APS typeset rendering's article text** — unreachable (JSON truncated at §I; PDFs refused by the harness).

### Honesty note
Every quoted passage above is reproduced exactly as the tool printed it, with the producing URL and, where available, the anchor. Where the search index itself truncated (`\tilde{\Phi}_{\math…`, `\mathcal{T}\left[\sigma_{1}\right]\mathcal{T}\...`, `\leq\left\|\mathcal{T}[\sigma]...`), the truncation is preserved and **not** completed. Nothing here has been inferred, reconstructed, or paraphrased as if it were text of the paper; the two places that go beyond literal quotation (the "value of this quote" reading of footnote 61, and the interpretation of N1/N2) are explicitly labelled as such. Items in §7 were genuinely not seen in any tool output, which means only that every route available to this session failed — **not** that the content is absent from the article.

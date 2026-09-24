# Section IV and Section V of arXiv:2212.11743v3 / Phys. Rev. X 14, 021011 — retrieval attempt

Target: Hong-Yi Wang, Fei Song, Zhong Wang, "Amoeba Formulation of Non-Bloch Band Theory in Arbitrary Dimensions",
arXiv:2212.11743v3, Phys. Rev. X 14, 021011 (2024), DOI 10.1103/PhysRevX.14.021011 (open access CC-BY; 21 pages).

Goal: verbatim text of Sec. IV ("Energy Spectra and density of states": IV.1 Statement of the proposal,
IV.2 Numerical evidence, IV.3 Derivation) and Sec. V ("Amoeba hole closing and spectral boundary").

**BOTTOM LINE (stated up front, no hedging): the running text of Sections IV and V could NOT be read.**
No route reachable from this session returned the body text of those sections. Everything below is either
(a) text I genuinely read from those sections' *neighbourhood* (Secs. I–III, plus the Sec. IV/V headings in the
table of contents), or (b) short verbatim *fragments* of the paper's later pages that were recovered indirectly
from web-search index snippets (each marked as such), or (c) clearly-labelled **other papers**, which are NOT
verbatim text of this paper. Nothing below is invented or paraphrased in place of a quote.

---

## 1. Route table (every URL actually attempted in this session)

| # | URL | Result | Which portion of the article it contained |
|---|-----|--------|-------------------------------------------|
| 1 | https://content.openalex.org/works/W4394855923.grobid-xml | HTTP 401 | none — "API key required" for content downloads |
| 2 | https://api.openalex.org/works/W4394855923?select=content_urls,has_content,best_oa_location,primary_location | HTTP 200 | metadata only (content_urls point back to keyed OpenAlex content; OA location = APS PDF) |
| 3 | https://arxiv.org/html/2212.11743v3/ | HTTP 200, **TRUNCATED** | Abstract, Sec. I, Sec. II (II.1, II.2) through Eq. (13); Sec. III only its first sentences. Cut mid-sentence at "The amoeba of f is defined as the log-moduli of the zero locus of f," |
| 4 | https://ar5iv.labs.arxiv.org/html/2212.11743 | HTTP 200, **TRUNCATED** | identical rendering and identical cut point as #3 |
| 5 | https://arxiv.org/html/2212.11743v3/main.html | HTTP 404 | none ("No HTML for '2212.11743v3'") |
| 6 | https://arxiv.org/pdf/2212.11743v3 | rejected | unsupported content type "application/pdf" — tool never returns PDF bodies |
| 7 | https://lite.duckduckgo.com/lite/?q=%22amoebic+spectrum%22 | network failure | none |
| 8 | https://scholar.archive.org/search?q=%22amoebic+spectrum%22 | network failure | none |
| 9 | https://www.mojeek.com/search?q=%22amoeba+hole+closing%22 | HTTP 403 | none (bot block) |
| 10 | https://www.bing.com/search?q=%22amoeba+hole+closing%22 | redirect to cn.bing.com, not followed | none |
| 11 | https://searx.be/search?q=%22amoeba+hole+closing%22&format=json | network failure | none |
| 12 | https://arxiv-org.ezproxy.obspm.fr/html/2212.11743v3 | redirect to login.ezproxy.obspm.fr | none |
| 13 | https://search.brave.com/search?q=%22amoebic+spectrum%22 | network failure | none |
| 14 | https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevX.14.021011/fulltext | rejected | returned **PDF**, unsupported content type |
| 15 | https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevX.14.021011 | HTTP 200 | full bibliographic metadata + verbatim abstract (Sec. I–IX absent) |
| 16 | https://core.ac.uk/search?q=%22Amoeba+Formulation+of+Non-Bloch+Band+Theory%22 | HTTP 403 | none (Cloudflare "Just a moment...") |
| 17 | https://www.themoonlight.io/en/review/amoeba-formulation-of-non-bloch-band-theory-in-arbitrary-dimensions | HTTP 429 | none (Vercel security checkpoint) |
| 18 | https://api.semanticscholar.org/graph/v1/snippet/search?query=amoebic%20spectrum%20Ronkin&limit=10 | HTTP 429 | none (rate limited, 4 attempts, all 429) |
| 19 | https://arxiv.org/html/2212.11743v1 | HTTP 200, **TRUNCATED** | v1 (titled "Amoeba formulation of the non-Hermitian skin effect in higher dimensions"): Sec. I, II, and Sec. III through the amoeba definition Eq. (9) and the opening of the Ronkin-function definition. Cut at "R_f(μ) = ∫_{T^d} (" |
| 20 | https://harvest.aps.org/v2/journals/articles/10.1103/5s1z-5t9r/fulltext | rejected | PDF (control test on a different article) |
| 21 | https://api.microlink.io/?url=…arxiv.org/html/2212.11743v3/…&data.sec.selector=%23S4 (CSS-selector slice of Sec. IV) | HTTP 429 | none — "Your daily rate limit has been reached" |
| 22 | https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevX.14.021011/fulltext?format=xml | rejected | still PDF |
| 23 | https://www.google.com/search?q=%22amoebic+spectrum%22+%22hole+closing%22&num=30 | network failure | none |
| 24 | https://www.arxiv-vanity.com/papers/2212.11743/ | network failure | none |
| 25 | https://r.jina.ai/https://arxiv.org/html/2212.11743v3/ | network failure | none |
| 26 | https://r.jina.ai/https%3A%2F%2Farxiv.org%2Fhtml%2F2212.11743v3%2F%23S4 (target-selector attempt) | network failure | none |
| 27 | https://www.emergentmind.com/papers/2212.11743 | HTTP 200, truncated | abstract only (rest behind sign-up) |
| 28 | https://api.fatcat.wiki/v0/release/lookup?doi=10.1103/physrevx.14.021011&expand=files | network failure | none |
| 29 | https://cn.bing.com/search?q=%22amoebic+spectrum%22&format=rss | HTTP 200 | RSS feed but region-biased, irrelevant results; no passage of the paper |
| 30 | https://www.textise.net/showText.aspx?strURL=…ar5iv…2212.11743 | HTTP 403 | none (Cloudflare block; egress IP reported as 122.97.209.4, China) |
| 31 | https://arxiv.org/html/2212.11743v3/index.html | HTTP 404 | none |
| 32 | https://www.instapaper.com/text?u=https%3A%2F%2Farxiv.org%2Fhtml%2F2212.11743v3%2F | HTTP 200 | login page only |
| 33 | https://xueshu.baidu.com/s?wd=%22amoebic%20spectrum%22 | HTTP 403 | none (百度安全验证 CAPTCHA) |
| 34 | https://www.bing.com/search?q=%22amoebic+spectrum%22&ensearch=1&format=rss | redirect to cn.bing.com | none |
| 35 | https://www.semanticscholar.org/arxiv/2212.11743 | HTTP 404 | none |
| 36 | https://papers.cool/arxiv/2212.11743 | HTTP 200 | abstract + metadata only |
| 37 | https://www.alphaxiv.org/abs/2212.11743 | HTTP 200, truncated | abstract + **third-party AI "Overview"** (see §4; not the paper's words) |
| 38 | https://www.x-mol.com/paper/1780375393846292480 | HTTP 200 | CAPTCHA page only |
| 39 | https://ouci.dntb.gov.ua/en/?s=10.1103%2FPhysRevX.14.021011 | HTTP 200 | search shell, no record rendered |
| 40 | https://www.alphaxiv.org/pdf/2212.11743 | HTTP 200 | JavaScript reader shell, no text |
| 41 | https://www.alphaxiv.org/docs/mcp | network failure | none |
| 42 | https://api.alphaxiv.org/v1/papers/2212.11743 | HTTP 404 | none |
| 43 | https://www.alphaxiv.org/api/paper/2212.11743 | timeout | none |
| 44 | https://arxiv.org/abs/2503.11505 | HTTP 200 | abstract of a **different** paper (Yang & Fang) |
| 45 | https://xxx.itp.ac.cn/html/2212.11743v3 | DNS ENOTFOUND | none |
| 46 | http://archive.org/wayback/available?url=arxiv.org/html/2212.11743v3 | network failure | none |
| 47 | https://journals.aps.org/prx/fulltext/10.1103/PhysRevX.14.021011/body | HTTP 404 | none — no separately addressable sub-document |
| 48 | https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevX.14.021011/fulltext?format=jats | rejected | still PDF |
| 49 | https://api.ocr.space/parse/imageurl?apikey=helloworld&url=…arxiv.org/pdf/2212.11743v3&page=9 | HTTP 400 | none — API rejects a "page" parameter; no per-page OCR available |
| 50 | https://www.w3.org/services/html2txt?url=https%3A%2F%2Farxiv.org%2Fhtml%2F2212.11743v3%2F | HTTP 429 | none (Cloudflare interstitial) |
| 51 | https://www.nature.com/articles/s42005-026-02546-2 | redirect to idp.nature.com (auth) | none |
| 52 | https://arxiv.org/html/2502.17931v1 | HTTP 200, truncated | **different paper** (Kaneshiro & Peters); Sec. II.3 "Amoeba formulation" read in full — see §4 |
| 53 | https://journals.aps.org/prx/supplemental/10.1103/PhysRevX.14.021011 | HTTP 404 | no supplemental material exists |
| 54 | https://arxiv.org/html/2212.11743v2 | HTTP 200, **TRUNCATED** | same cut point as v3 (end of Sec. II.2, Eq. (13)) |
| 55 | https://api.crossref.org/works/10.1103/PhysRevX.14.021011 | HTTP 200 | metadata, reference list (56 refs), abstract |
| 56 | https://journals.aps.org/prx/xml/10.1103/PhysRevX.14.021011 | network failure | none |
| 57 | https://arxiv.org/html/2407.01296v1 | HTTP 200, truncated | **different paper** (Xiong, Xing & Hu); Sec. I–II.2 read — see §4 |
| 58 | https://html.duckduckgo.com/html/?q=%22amoebic+spectrum%22+%22hole+closing%22 | network failure | none |
| 59 | https://search.marginalia.nu/search?query=%22amoebic+spectrum%22 | redirect → marginalia-search.com | none |
| 60 | https://www.startpage.com/sp/search?query=%22amoebic+spectrum%22 | network failure | none |
| 61 | https://browse.arxiv.org/html/2212.11743v3 | DNS ENOTFOUND | none |
| 62 | https://marginalia-search.com/search?query=%22amoebic+spectrum%22 | HTTP 200 | bot-throttle page + help text, no results |
| 63 | https://www.ecosia.org/search?q=%22amoebic+spectrum%22+%22hole+closing%22 | redirect to www.bing.com | none |

Search-side queries issued (the built-in search tool's index DOES hold the publisher PDF full text and returned
verbatim passage fragments; these queries are the source of every fragment in §3). Recorded for reproducibility:
"amoebic spectrum"; "amoeba hole closing"; "amoebic spectrum" "Ronkin function"; "the minimum of the Ronkin function";
"we propose" "amoebic spectrum"; "disk geometry" "square geometry"; "random geometry" amoeba;
"universal spectrum and universal DOS"; "the proof of this theorem is now simple"; "Taking a circle" "that surrounds
the whole PBC spectrum"; "The language of Toeplitz matrices is very useful in addressing tight-binding Hamiltonians";
"the precise meaning of which will be discussed"; "amoebic spectrum" "Fig."; "amoebic spectrum" "theorem";
"hole closes" amoeba; "spectral boundary" amoeba; "Statement of the proposal" amoeba; "Numerical evidence" amoeba;
"the spectral potential" amoeba "density of states"; plus ~50 more variants. All returned the same small set of
passages (reproduced in §3); none returned the body of Sec. IV.1 or IV.2.

---

## 2. What I actually read from the article itself (verbatim)

### 2.1 Section headings of the target sections, as they appear in the article (verbatim, from the arXiv HTML table of contents)

```
4.  [IV Energy Spectra and density of states]
    1.  [IV.1 Statement of the proposal]
    2.  [IV.2 Numerical evidence]
    3.  [IV.3 Derivation]
5.  [V Amoeba hole closing and spectral boundary]
```

Also from the same rendered article (v3), the abstract and the two introduction sentences that forward-reference
Sections IV and V verbatim:

> The non-Hermitian skin effect dramatically reshapes the energy bands of non-Hermitian systems, meaning that the
> usual Bloch band theory is fundamentally inadequate as their characterization. The non-Bloch band theory, in which
> the concept of Brillouin zone is generalized, has been widely applied to investigate non-Hermitian systems in one
> spatial dimension. However, its generalization to higher dimensions has been challenging. Here, we develop a
> formulation of the non-Hermitian skin effect and non-Bloch band theory in arbitrary spatial dimensions, which is
> based on a natural geometrical object known as the amoeba. Our theory provides a general framework for studying
> non-Hermitian bands beyond one dimension. Key quantities of non-Hermitian bands, including the energy spectrum,
> eigenstates profiles, and the generalized Brillouin zone, can be efficiently obtained from this approach.

> We show in a theorem that the energy spectrum can be obtained from the shape of the amoeba. We also demonstrate,
> despite the geometry-dependent NHSE, the existence of a universal spectrum (amoebic spectrum) to which the OBC
> spectrum under any generic geometry converges.

> In Secs. [IV]–[VI], we introduce the amoeba formulation for non-Hermitian systems, and then make use of this
> formulation and the theory of Toeplitz matrices to establish a universal way to determine the DOS as well as the
> GBZ.

The same forward reference in the v1 preprint reads (v1 title: "Amoeba formulation of the non-Hermitian skin effect
in higher dimensions"):

> A rigorous theorem is proved as the general basis for DOS calculations. We also demonstrate, despite the
> geometry-dependent NHSE, the existence of a universal spectrum to which the OBC spectrum under a generic geometry
> converges.

### 2.2 The last text before the truncation (Sec. II.2 end — the observation that Sec. IV builds on) — verbatim

> From the above examples, we observe that the absence (presence) of a hole in the amoeba of the characteristic
> polynomial could be an indicator of the energy E being (not being) in the OBC energy spectrum. This is a key
> observation of the present work. To obtain more quantitative results from this observation, it is helpful to know
> some mathematical properties about the amoeba.

### 2.3 Sec. III opening (all that is visible of Sec. III in v3) — verbatim, then TRUNCATED

> In this section, we shall introduce the basic concept of amoeba and a closely related analytic tool, the Ronkin
> function. As a quite recent concept in mathematics, the amoeba was introduced by Gelfand et al. in 1994. Albeit
> elementary, the notion of amoeba has deep connections with various concepts in algebraic geometry, which has
> stimulated extensive studies in mathematics.
>
> Let f be a Laurent polynomial of β_j, j=1,2,…,d, where d will be identified as the spatial dimension in our study.
> The amoeba of f is defined as the log-moduli of the zero locus of f,
>
> *(the tool's output was cut off here: "The amoeba of f is defined as the log-moduli of the zero locus of f,")*

From the **v1** HTML (same section, slightly different numbering) I additionally read, verbatim:

> The amoeba of f is defined as the log-moduli of the zero locus of f
>
> A_f = { log|β| : f(β)=0 } ⊂ ℝ^d,   (9)

and the first line of the Ronkin-function definition:

> A useful analytic tool in the study of amoeba is the Ronkin function, which is defined as [27, 31, 32]:
>
> R_f(μ) = ∫_{T^d} (   *(cut off here)*

Nothing further in either version could be retrieved. **The Ronkin function's explicit formula as written by these
authors (the definition on which the Sec. IV potential φ(E) is built) therefore could not be read.**

---

## 3. Verbatim fragments of the LATER pages (Secs. IV.3 / V region), recovered only as search-index snippets

These are the only pieces of the paper's later text I obtained. They come from the search tool's index of the
publisher PDF (and of the arXiv/ar5iv HTML). Every fragment below is reproduced **exactly as the tool returned it**;
the trailing `...` is the *tool's* truncation, not mine. I do **not** know the equation numbers, and I cannot place
these fragments at exact line/section positions. Fragments 1–6 are the recurring set; the anchors shown were
attached by the search tool to the source PDF URL.

1. (source: https://journals.aps.org/prx/pdf/10.1103/PhysRevX.14.021011, snippet anchor `#9#3`)
   > We now apply the explicit formula of the Ronkin function to \(\operatorname *{det}[E - h(\beta)] = a_{- M}(E)\beta^{- M} + \dots...

2. (source: same PDF, snippet anchor `#9#4`; also returned against the ar5iv HTML and the ezproxy HTML mirror)
   > We emphasize that the Ronkin function tells the unique universal spectrum and universal DOS of the OBC system, the precise meani...

   A variant of this same sentence, returned against the ar5iv copy (no "unique"), reads:
   > We emphasize that the Ronkin function tells the universal spectrum and universal DOS of the OBC system, the precise meaning of w...

3. (source: same PDF, snippet anchor `#9#4`)
   > The language of Toeplitz matrices is very useful in addressing tight- binding Hamiltonians

4. (source: same PDF, snippet anchor `#9#5`)
   > With the results of the DOS in the previous section, the proof of this theorem is now simple

5. (source: same PDF, snippet anchor `#9#7`)
   > Taking a circle \(|E| = R\) that surrounds the whole PBC spectrum, the average potential on this circle is \(\tilde{\Phi}_{\math...

6. (source: same PDF, snippet anchor `#9#8`)
   > Therefore, a condition for Eq

7. (source: http://export-test.arxiv.org/pdf/2212.11743, snippet anchor `#9#7`)
   > An alternative proof of the spectral inequality is based on Eq

8. (source: http://export-test.arxiv.org/pdf/2212.11743, snippet anchor `#9#3`)
   > \[R_{g}(\log|R|) =\log|g(0)|+\log\frac{|z_{2}|}{|z_{1}|}+2\log\frac{|z_{3}|}{|z_{2}| }+3\log\frac{|z_{4}|}{|z_{3}|}\] \[+\cdots+...

9. (source: the ar5iv/ezproxy HTML of 2212.11743v3, snippet anchor `#7`)
   > For the present model, the precise phase boundary can even be analytically determined by the amoeba formulation

10. (source: same PDF, snippet anchor `#9#6`)
   > This model is known to have a Chern insulator phase as well as a trivial insulator phase in the Hermitian limit \((\gamma = 0)\)...

11. (source: same PDF, snippet anchor `#9#2`) — the model Hamiltonian, i.e. Eq. (12) of Sec. II.2, recovered as a
    later-page fragment:
    > \[\begin{array}{c}{h(\beta) = t(\beta_x + \beta_x^{-1} + \beta_y + \beta_y^{-1})}\\ {+t'(\beta_x + \beta_x^{-1})(\beta_y + \beta...

**These fragments are NOT a statement of the proposal.** In particular, I could NOT obtain:
the wording of any "Theorem", the definition of φ(E) as the minimum of the Ronkin function, the DOS formula
ρ(E) = (1/2π)Δφ(E) as printed in this paper, the "amoebic spectrum" definition statement, any equation number from
Sec. IV (or from Sec. V), any statement in Sec. IV.2 about disk/square/random geometry or system sizes, the
derivation outline of Sec. IV.3, or the Sec. V statements about hole closing / the spectral boundary.
Searching the exact phrase "the minimum of the Ronkin function" returned **no hit inside this paper**.

---

## 4. Secondary sources — explicitly NOT verbatim text of this paper

These are **other documents**. They are included only because they restate the same formalism, and they are labelled
as such. Do not attribute any of the following to Wang–Song–Wang except where a source says so.

### 4.1 alphaXiv "AI Overview" for arXiv:2212.11743 (machine-generated third-party summary, NOT the paper's words)
Source: https://www.alphaxiv.org/abs/2212.11743 (truncated). It states (as its own generated prose):
- the amoeba A_f = {(log|β_1|, …, log|β_d|) ∈ ℝ^d | f(β_1,…,β_d) = 0};
- the Ronkin function R_f(μ) = (1/(2π)^d) ∫_{[0,2π]^d} log|f(e^{μ_1+iθ_1},…,e^{μ_d+iθ_d})| dθ_1…dθ_d;
- "ϕ(E) = min_μ R_{det(E−h)}(μ)";
- "ρ(E) = (1/2π) Δϕ(E) where Δ = ∂²/∂(ReE)² + ∂²/∂(ImE)²";
- "The OBC spectrum corresponds exactly to the regions in the complex energy plane where this central hole in the
  amoeba has closed."
This is AI-generated text about the paper; I record it as corroboration of *content*, never as a quotation of the paper.

### 4.2 Kaneshiro & Peters, "Symplectic-Amoeba formulation of the non-Bloch band theory for one-dimensional two-band
systems", arXiv:2502.17931, Sec. II.3 "Amoeba formulation" (read verbatim from
https://arxiv.org/html/2502.17931v1). Their equations, verbatim as rendered:

> ρ(E) = lim_{N→∞} (1/N) Tr δ(E − H^{(N)}).   (7)
>
> ρ(E) = (1/2π) Δ ϕ(E),   where Δ = ∂²/∂(Re E)² + ∂²/∂(Im E)².   (8)
>
> ϕ(E) = lim_{N→∞} (1/N) Tr ln|E − H^{(N)}|.   (9)
>
> ϕ(E) = ∫_{ln|β|=0} (dβ / 2πiβ) ln|ChP(β, E)|,   (10)   [one-band class A, Szegő's limit theorem]
>
> Φ(E) = min_μ R_E(μ),   (12)
>
> R_E(μ) = ∮_{ln|β|=μ} (dβ / 2πiβ) ln|ChP(E, β)|.   (13)
>
> R_E(μ) = ln|C_E| − pμ + Σ_{j=1}^{2p} μ_j + Σ_{j=1}^{2p} (μ − μ_j)·θ(μ − μ_j),   (15)
>
> "In the Amoeba formulation, the absence of an intermediate region where the Ronkin function is constant indicates
> that E lies within the spectrum, and the μ*, which minimizes the function, yields the corresponding eigenstate's
> inverse localization length. This criterion, known as "hole closing," is consistent with the GBZ condition..."

(Their Ref. [58] is the Wang–Song–Wang paper, i.e. this is their restatement, not the original wording.)

### 4.3 Xiong, Xing & Hu, arXiv:2407.01296, Sec. II.2 (read verbatim from https://arxiv.org/html/2407.01296v1):

> ρ(E) = (1/N) Σ_{i=1}^{N} δ_{E,E_i},   (9)
>
> ϕ(E) = (1/N) Σ_{i=1}^{N} log|E_i − E|,   (10)
>
> ρ(E) = (1/2π) ∇²_E ϕ(E),   (11)
>
> lim_{N→∞} ϕ(E) = ∫_0^{2π} (dk/2π) log|det[H(e^{ik+μ}) − E]|,   (12)

Their Sec. V.1 is titled "Amoeba formulation" and their text says of the amoeba formulation that it "neglects
geometric information and yields geometry-irrelevant non-Bloch spectra" (Ref. [65] there = the Wang–Song–Wang paper).

### 4.4 Supplemental Material of the 2026 paper "Geometry-adaptive formulation of non-Bloch bands in arbitrary
dimensions and spectral instability" (search snippet only, source
https://static-content.springer.com/esm/art%3A10.1038%2Fs42005-026-02546-2/MediaObjects/42005_2026_2546_MOESM2_ESM.pdf):

> \[\begin{array}{rl} & {\Phi_{\mathrm{Amoeba}}(E) = \iint \frac{dk_1dk_2}{(2\pi)^2}\log |f(e^{ik_1 + \mu_1,\min},e^{ik_2 + \mu_2,...

### 4.5 Abstract as served by the APS Harvest API and Crossref (verbatim, metadata service):
https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevX.14.021011 —
authors "Hong-Yi Wang", "Fei Song", "Zhong Wang"; affiliation "Institute for Advanced Study, Tsinghua University,
Beijing 100084, China"; journal PRX vol. 14, issue 2, pageStart 021011, numPages 21, date 2024-04-16, CC-BY 4.0,
arXiv:2212.11743. Abstract text identical to §2.1 above.

---

## 5. Everything I could NOT read (explicit list)

1. **Sec. IV.1 "Statement of the proposal"** — 100% not read. No statement of the proposal, no theorem wording, no
   φ(E) = min Ronkin definition, no ρ(E) = (1/2π)Δφ(E) formula, no "amoebic spectrum" definition, no equation numbers.
2. **Sec. IV.2 "Numerical evidence"** — 100% not read. No statements about disk / square / random geometry, no system
   sizes, no convergence figures or captions (Figs. 3–5 region not read).
3. **Sec. IV.3 "Derivation"** — not read. Only the isolated fragment "We now apply the explicit formula of the Ronkin
   function to det[E − h(β)] = a_{−M}(E)β^{−M} + …" and "The language of Toeplitz matrices is very useful in
   addressing tight-binding Hamiltonians" (Sec. IV.3 or V region), plus "Therefore, a condition for Eq".
4. **Sec. V "Amoeba hole closing and spectral boundary"** — not read. Only the isolated fragments in §3
   (items 4, 5, 6, 8, 9). No hole-closing criterion, no spectral-boundary characterisation, no equation numbers.
5. **Sec. III's mathematical content beyond its first paragraph** — not read (truncated at the amoeba definition;
   the Ronkin-function formula as printed by these authors was cut off mid-expression in v1 as well).
6. **Any figure, table, caption, or appendix of this paper** — not read, except the Fig. 1/Fig. 2 captions that fall
   in Sec. I–II (which I did read; not reproduced here since they are outside the target sections).
7. **The publisher's HTML/JSON full text (journals.aps.org/prx/fulltext/10.1103/PhysRevX.14.021011)** — the head-only
   truncation limit means it yields Sections I–III only; a direct fetch was therefore not useful, and its
   `/body` and `/supplemental` sub-paths return 404.
8. **APS Harvest full text** — exists but is served as `application/pdf`, which this tool refuses; no JATS/XML variant
   responded (`?format=xml`, `?format=jats` both return PDF).
9. **OpenAlex GROBID XML / PDF content** — HTTP 401, requires an API key not available in this session.
10. **Semantic Scholar snippet-search API** — HTTP 429 on all four attempts (no key).
11. **microlink CSS-selector slicing of `#S4`** — HTTP 429 (daily limit reached); this was the one route that would
    have returned exactly Section IV.
12. **PDF text extraction of any kind** — no route worked: direct PDF fetch is refused (content type), OCR.space has
    no per-page parameter, the W3C html2txt service is Cloudflare-blocked, and no keyless PDF→text API was reachable.
13. **Google/Bing/Mojeek/Startpage/Ecosia/Brave/Baidu/DuckDuckGo/Marginalia result pages** — blocked, redirected, or
    useless, due to the session's Chinese egress IP (Cloudflare report: 122.97.209.4).
14. **r.jina.ai, arxiv-vanity, scholar.archive.org, CORE, themoonlight.io, Instapaper text view, archive.org/fatcat,
    xxx.itp.ac.cn** — network failures, auth walls, or 403/429.

**Honest summary:** with the tools available in this session the later half of this paper is not retrievable.
The full-body routes (OpenAlex keyed content, APS Harvest XML, r.jina.ai, text proxies, PDF parsers, CSS-selector
extraction) are each blocked by a distinct, verified obstacle, and the only leakage of late-section text is the set of
~11 search-index fragments quoted verbatim in §3.

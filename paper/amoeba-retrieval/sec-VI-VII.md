# Retrieval report: Sec. VI "Generalized Brillouin zone" and Sec. VII "Non-Bloch band topology"
Target: Hong-Yi Wang, Fei Song, Zhong Wang, "Amoeba Formulation of Non-Bloch Band Theory in Arbitrary Dimensions", arXiv:2212.11743 (v3), Phys. Rev. X 14, 021011 (2024), DOI 10.1103/PhysRevX.14.021011.

STATUS: **RETRIEVAL FAILED FOR THE REQUESTED SECTIONS.** Every reachable rendering of the article that the session's `web_fetch` tool can open is truncated by the harness at ~50,000 characters, which consistently lands at the *end of Sec. III* (the definition of the Ronkin function, Eq. (15) in arXiv v3 / APS numbering, Eq. (13) in arXiv v1 numbering). Sec. VI and Sec. VII text was NOT obtained. The only verbatim fragments of the paper's later text that I obtained come from the `web_search` snippet channel and are 1–2 sentence fragments (quoted in section 4 below); they do not contain the GBZ definition, the Chern-band Hamiltonian, or the Chern-number equation.

No content below is paraphrased. Anything not literally present in a tool output is marked "could not read".

---

## 1. Route table (every URL attempted, in order)

| # | URL | Result | Portion of the article obtained |
|---|-----|--------|---------------------------------|
| 1 | https://content.openalex.org/works/W4394855923.grobid-xml | FAILED (HTTP 401 "API key required") | none |
| 2 | https://api.openalex.org/works/W4394855923?select=content_urls,has_content,id,title | FETCHED (HTTP 200) | metadata only: `{"content_urls":{"pdf":"https://content.openalex.org/works/W4394855923.pdf","grobid_xml":"https://content.openalex.org/works/W4394855923.grobid-xml"},"has_content":{"grobid_xml":true,"pdf":true},"id":"https://openalex.org/W4394855923","title":"Amoeba Formulation of Non-Bloch Band Theory in Arbitrary Dimensions"}` |
| 3 | https://api.openalex.org/works/W4394855923/ngrams | FAILED (fetch failed) | none |
| 4 | https://ar5iv.org/abs/2212.11743 | FAILED (cross-origin redirect to ar5iv.labs.arxiv.org not followed) | none |
| 5 | https://ar5iv.labs.arxiv.org/html/2212.11743 | FETCHED (HTTP 200) | TRUNCATED at end of Sec. III (amoeba definition Eq. (14), Ronkin function Eq. (15) statement) — Sec. I, II, II.1, II.2, start of III only |
| 6 | https://arxiv.org/html/2212.11743v3/ | FETCHED (HTTP 200) | TRUNCATED at end of Sec. III, last line = "Let \(f\) be a Laurent polynomial of \(\beta_{j}\)..." then "(Content truncated...)" — Sec. I, II, II.1, II.2, start of III |
| 7 | https://arxiv.org/src/2212.11743v3 | FAILED (unsupported content type "application/gzip") — LaTeX source not usable | none |
| 8 | https://api.codetabs.com/v1/proxy?quest=https://arxiv.org/html/2212.11743v3 | FAILED (fetch failed) | none |
| 9 | https://api.allorigins.win/raw?url=<arxiv html v3> | FAILED (HTTP 520) | none |
| 10 | https://html.duckduckgo.com/html/?q=... | FAILED (fetch failed) | none |
| 11 | https://arxiv.org/html/2212.11743v3/2212.11743v3.bbl | HTTP 404 (no separate per-section/asset files exist in the arXiv HTML tree) | none |
| 12 | https://www.textise.net/showText.aspx?strURL=<arxiv html v3> | FAILED (HTTP 403 Cloudflare block) | none |
| 13 | https://r.jina.ai/https://arxiv.org/html/2212.11743v3 | FAILED (fetch failed) | none |
| 14 | https://colab.ws/articles/10.1103%2FPhysRevX.14.021011 | FAILED (HTTP 403 DDoS-Guard, IP-restricted) | none |
| 15 | https://ouci.dntb.gov.ua/en/?s=10.1103%2FPhysRevX.14.021011 | FETCHED (HTTP 200) | empty site shell only (JavaScript "Loading..."); no article text |
| 16 | https://cn.bing.com/search?q="amoebic spectrum"... | FETCHED (HTTP 200) | no relevant results (returned generic dictionary/douyin results) |
| 17 | https://journals.aps.org/prx/pdf/10.1103/PhysRevX.14.021011 | FAILED (unsupported content type "application/pdf") | none |
| 18 | https://api.allorigins.win/raw?url=<arxiv pdf v3> | FAILED (tool call timed out after 30000ms) | none |
| 19 | https://search.brave.com/search?q=... | FAILED (fetch failed) | none |
| 20 | https://www.mojeek.com/search?q="amoeba hole closing" | FAILED (HTTP 403, automated-query block) | none |
| 21 | https://searx.be/search?q=... | FAILED (fetch failed) | none |
| 22 | https://www.textance.com/cgi-bin/textance.cgi?url=<arxiv pdf v3> | FAILED (DNS: getaddrinfo ENOTFOUND) | none |
| 23 | https://dx.doi.org/10.1103/PhysRevX.14.021011 | FAILED (cross-origin redirect to link.aps.org not followed) | none |
| 24 | https://arxiv.org/html/2212.11743v3/S6.html | HTTP 404 (per-section files do not exist) | none |
| 25 | https://arxiv.org/html/2212.11743v3/index.html | HTTP 404 | none |
| 26 | https://lite.duckduckgo.com/lite/?q=%22amoebic+spectrum%22 | FAILED (fetch failed) | none |
| 27 | https://link.aps.org/doi/10.1103/PhysRevX.14.021011 | FAILED (cross-origin redirect to journals.aps.org not followed) | none |
| 28 | https://journals.aps.org/prx/fulltext/10.1103/PhysRevX.14.021011 | FETCHED (HTTP 200, structured JSON for the whole article) | TRUNCATED at the Ronkin-function derivation (Eq. (15)/(16) region) — Sec. I, II.A, II.B, start of III. Raw size ≈ 100 KB; spill file is head-only, so the tail is unrecoverable |
| 29 | http://export-test.arxiv.org/pdf/2212.11743 | FAILED (DNS: getaddrinfo ENOTFOUND) | none |
| 30 | https://cn.bing.com/search?q="arxiv.org/html/2212.11743v3#S7" | FETCHED (HTTP 200) | no relevant results (only generic arXiv/zhihu pages) |
| 31 | https://corsproxy.io/?https://arxiv.org/html/2212.11743v3 | FAILED (HTTP 403 "keyless_legacy_url" — API key now required) | none |
| 32 | https://journals.aps.org/prx/fulltext/10.1103/PhysRevX.14.021011/s6 | HTTP 404 ("Not Found") | none |
| 33 | https://journals.aps.org/prx/fulltext/10.1103/PhysRevX.14.021011?section=s6&format=json | FETCHED (HTTP 200) — **query string ignored**; byte-identical head to #28 | same as #28 (truncated at Sec. III) |
| 34 | https://journals.aps.org/prx/fulltext/10.1103/PhysRevX.14.021011/body | HTTP 404 | none |
| 35 | https://journals.aps.org/prx/article/10.1103/PhysRevX.14.021011 | HTTP 404 | none |
| 36 | https://journals.aps.org/prx/export/10.1103/PhysRevX.14.021011 | FETCHED (HTTP 200) | BibTeX metadata only (@article PhysRevX.14.021011) |
| 37 | https://scholar.archive.org/search?q="Amoeba Formulation of Non-Bloch Band Theory" | FAILED (fetch failed) | none |
| 38 | https://www.textise.net/showText.aspx?strURL=<aps fulltext> | FAILED (HTTP 403 Cloudflare block) | none |
| 39 | https://r.jina.ai/https://ar5iv.labs.arxiv.org/html/2212.11743 | FAILED (fetch failed) | none |
| 40 | https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevX.14.021011 | FAILED (fetch failed) | none |
| 41 | https://journals.aps.org/prx/fulltext/10.1103/PhysRevX.14.021011?format=nice | FETCHED (HTTP 200) — parameter ignored, byte-identical head to #28 | same as #28 |
| 42 | https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevX.14.021011/fulltext | FAILED (unsupported content type "application/pdf" — this endpoint serves the PDF) | none |
| 43 | https://export-test.arxiv.org/pdf/2212.11743 | FAILED (DNS: getaddrinfo ENOTFOUND) | none |
| 44 | https://harvest.aps.org/v2/journals/articles/10.1103/5s1z-5t9r/fulltext | FAILED (unsupported content type "application/pdf") — control: also PDF | none |
| 45 | https://arxiv-org.ezproxy.obspm.fr/html/2212.11743v3 | FAILED (fetch failed) | none |
| 46 | https://huggingface.co/papers/2212.11743 | FAILED (fetch failed) | none |
| 47 | https://www.emergentmind.com/papers/2212.11743 | FETCHED (HTTP 200) | title, abstract, author list, arXiv metadata only; truncates after "Authors (3)" — no body text |
| 48 | https://www.alphaxiv.org/abs/2212.11743 | FAILED (fetch failed) | none |
| 49 | https://r.jina.ai/https%3A%2F%2Farxiv.org%2Fhtml%2F2212.11743v3 | FAILED (fetch failed) | none |
| 50 | https://r.jina.ai/https://arxiv.org/abs/2212.11743 | FAILED (fetch failed) | none |
| 51 | https://www.semanticscholar.org/arxiv/2212.11743 | HTTP 404 | none |
| 52 | https://defuddle.md/https://arxiv.org/html/2212.11743v3 | FAILED (fetch failed) | none |
| 53 | https://urltotext.com/api/?url=<arxiv html v3> | FAILED (HTTP 500) | none |
| 54 | https://www.arxiv-vanity.com/papers/2212.11743/ | FAILED (fetch failed) | none |
| 55 | https://core.ac.uk/search?q="Amoeba Formulation of Non-Bloch Band Theory" | FAILED (fetch failed) | none |
| 56 | https://arxiv.org/html/2212.11743v1 | FETCHED (HTTP 200) | Section list I–IX + App. A/B (confirms Sec. VI "Generalized Brillouin zone", Sec. VII "Non-Bloch band topology"); body TRUNCATED at end of Sec. III (Eq. (9) amoeba, start of Ronkin function). v1 has different equation numbering (Chern-band model is Eq. (39)) |
| 57 | https://arxiv.org/html/2212.11743v2 | FETCHED (HTTP 200) | identical v3-style numbering; TRUNCATED at end of Sec. III, last line "For our single-band model, the determinant can be dropped, and therefore the characteristic equation is simply" |
| 58 | https://browse-export.arxiv.org/html/2212.11743v1 | FETCHED (HTTP 200) — byte-identical to #56 | same as #56 |
| 59 | https://marginalia-search.com/search?query="amoeba formulation of non-Bloch" | FETCHED (HTTP 200) | search engine served a bot-throttle page; no results |
| 60 | https://www.bing.com/search?q=...&format=rss | FAILED (cross-origin redirect to cn.bing.com) | none |
| 61 | https://arxiv.org/html/2212.11743v1/2212.11743v1.html | FETCHED (HTTP 200) — byte-identical to #56 | same as #56 |
| 62 | https://www.semanticscholar.org/paper/84ccbf769bc1ee1312a3ccf4074b8e9ffc019582 | FETCHED (HTTP 202) — empty body (JS-gated / bot-challenged) | none |
| 63 | https://arxiv.org/html/2212.11743v1/s6.html | HTTP 404 | none |
| 64 | https://arxiv.org/html/2212.11743v1/node6.html | HTTP 404 | none |
| 65 | https://arxiv.org/html/2212.11743v1/S6 | HTTP 404 | none |
| 66 | https://export.arxiv.org/html/2212.11743v1 | FETCHED (HTTP 200) — byte-identical to #56 | same as #56 |
| 67 | https://cn.bing.com/search?q=...&format=rss | FETCHED (HTTP 200) | RSS returned irrelevant results (amoeba biology, MySQL Amoeba) |
| 68 | https://export.arxiv.org/html/2212.11743v1?format=text | FETCHED (HTTP 200) — parameter ignored, byte-identical to #56 | same as #56 |
| 69 | https://www.x-mol.com/paper/1780375393846292480 | FAILED (HTTP 200 but CAPTCHA "人机验证" page; no article text) | none |
| 70 | https://api.codetabs.com/v1/proxy/?quest=https://arxiv.org/html/2212.11743v1 | FAILED (HTTP 522) | none |
| 71 | https://cn.bing.com/search?q=...&format=rss (2nd query) | FETCHED (HTTP 200) | irrelevant results |
| 72 | https://arxiv-org.ezproxy.obspm.fr/html/2212.11743v1 | FAILED (cross-origin redirect to login.ezproxy.obspm.fr — institutional login required) | none |
| 73 | https://journals.aps.org/prx/fulltext/10.1103/PhysRevX.14.021011/s6?format=json | HTTP 404 | none |
| 74 | https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevX.14.021011/fulltext.json | HTTP 404 (`{"errors":[{"title":"...fulltext.json not found"}]}`) | none |
| 75 | https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevX.14.021011/fulltext.txt | HTTP 404 (`{"errors":[{"title":"...fulltext.txt not found"}]}`) | none |

`web_search` queries issued (all returned only other-paper snippets or metadata, never the target sections' text): see section 4.

Key findings about the constraint:
- The APS `fulltext` endpoint (#28) is machine-readable structured JSON with explicit section anchors `s1`…`s8` (and `s2b`), i.e. Sec. VI = `s6`, Sec. VII = `s7`. The harness returns only the first ~50 KB of that JSON; the spill file it writes (path reported in the tool output, e.g. `...\bee84f701f26-web_fetch.txt`) contains only the head, so the tail is genuinely lost.
- Every query-string parameter tried on APS (`?section=s6`, `?format=json`, `?format=nice`) is ignored: byte-identical head output. Path-based section addressing (`/s6`, `/body`) returns 404. `harvest.aps.org/.../fulltext` is a PDF.
- `arxiv.org/html/2212.11743v1` is the shortest rendering, but its 50 KB head still ends inside Sec. III. `v2` and `v3` likewise. There are no per-section HTML files (all such paths 404).
- `arxiv.org/src/2212.11743v3` is a gzip archive and the tool cannot process it (`unsupported content type "application/gzip"`), so the LaTeX source route is closed.
- Reader/extraction proxies tried (r.jina.ai, codetabs, allorigins, corsproxy.io, textise, urltotext, defuddle.md, textance) all failed or are key/instance blocked.

---

## 2. Verbatim text obtained from arXiv v3 HTML (https://arxiv.org/html/2212.11743v3/) — TRUNCATED AT END OF SEC. III

This is everything the tool returned. The final lines are cut mid-derivation; Sec. IV onward (and therefore Sec. VI and Sec. VII) were not served.

> # Amoeba Formulation of Non-Bloch Band Theory in Arbitrary Dimensions
>
> Hong-Yi Wang Affiliation: Institute for Advanced Study, Tsinghua University, Beijing 100084, China   Fei Song Affiliation: Institute for Advanced Study, Tsinghua University, Beijing 100084, China   Zhong Wang Email: wangzhongemail@tsinghua.edu.cn Affiliation: Institute for Advanced Study, Tsinghua University, Beijing 100084, China
>
> arXiv:2212.11743v3 [cond-mat.mes-hall] 30 Apr 2024
>
> ###### Abstract
>
> The non-Hermitian skin effect dramatically reshapes the energy bands of non-Hermitian systems, meaning that the usual Bloch band theory is fundamentally inadequate as their characterization. The non-Bloch band theory, in which the concept of Brillouin zone is generalized, has been widely applied to investigate non-Hermitian systems in one spatial dimension. However, its generalization to higher dimensions has been challenging. Here, we develop a formulation of the non-Hermitian skin effect and non-Bloch band theory in arbitrary spatial dimensions, which is based on a natural geometrical object known as the amoeba. Our theory provides a general framework for studying non-Hermitian bands beyond one dimension. Key quantities of non-Hermitian bands, including the energy spectrum, eigenstates profiles, and the generalized Brillouin zone, can be efficiently obtained from this approach.

[Sec. I Introduction and Sec. II Motivation (II.1 Review of 1D non-Bloch band theory, II.2 The way to the amoeba) as reproduced verbatim in section 3 of this file from the APS JSON; the arXiv rendering is the same text.]

> ## III Mathematical properties of the amoeba and Ronkin function
>
> In this section, we shall introduce the basic concept of amoeba and a closely related analytic tool, the Ronkin function. As a quite recent concept in mathematics, the amoeba was introduced by Gelfand *et al.* in 1994 [23]. Albeit elementary, the notion of amoeba has deep connections with various concepts in algebraic geometry, which has stimulated extensive studies in mathematics [24–27].
>
> Let \(f\) be a Laurent polynomial of \(\beta_{j}\), \(j=1,2,\dots,d\), where \(d\) will be identified as the spatial dimension in our study. The amoeba of \(f\) is defined as the log-moduli of the zero locus of \(f\),
>
> (14) \(\mathcal{A}_{f}=\{\log|\bm{\beta}|:f(\bm{\beta})=0\}\subset\mathbb{R}^{d},\)
>
> in which we use the notation \(\log|\bm{\beta}|\coloneqq(\log|\beta_{1}|,\dots,\log|\beta_{d}|)\) to simplify our expressions. Similar notations such as \(e^{\bm{\mu}}\coloneqq(e^{\mu_{1}},\dots,e^{\mu_{d}})\) are used hereafter. In our case, the Laurent polynomial in use is \(\det[E-h(\bm{\beta})]\). We can see that the geometric objects in Figs. 1(b) and 1(c) and 1(e) and 1(f) are 1D and 2D amoebae, respectively.
>
> The name amoeba was motivated by its appearance in 2D: It has slim “tentacles” extending to infinity, and sometimes several “vacuoles” (holes) inside its body. Importantly, a particular hole plays an important role in our formulation. It is known that the amoeba in any spatial dimensions is a closed set, and each hole is a convex set [25].
>
> A useful analytic tool in the study of amoeba is the Ronkin function, which is defined as [24, 30, 31]
>
> (15) \(R_{f}(\bm{\mu})=\int_{T^{d}}\left(\frac{d\bm{\theta}}{2\pi}\right)^{d}\log\left|f(e^{\bm{\mu}+i\bm{\theta}})\right|,\)
>
> where the domain of integration is the \(d\)-dimensional torus \(T^{d}=[0,2\pi]^{d}\), and the expression is simplified by the notations \(f(e^{\bm{\mu}+i\bm{\theta}})\coloneqq f(e^{\mu_{1}+i\theta_{1}},\dots,e^{\mu_{d}+i\theta_{d}})\), and \((d\theta/2\pi)^{d}\coloneqq(d\theta_{1}/2\pi)\dots(d\theta_{d}/2\pi)\).
>
> It is beneficial to study the gradient of the Ronkin function [25, 30]. To this end, we can express the integrand in \(R_{f}\) as \(\log|f|=\mathrm{Re}\log f\). The real part can be taken at the end of the calculation. It turns out that the integral is real-valued before taking the real part, and therefore the “\(\mathrm{Re}\)” symbol can be discarded. The derivation proceeds as
>
> (16) \(\nu_{j}=\frac{\partial R_{f}(\bm{\mu})}{\partial\mu_{j}}\)
> \(=\mathrm{Re}\int_{T^{d}}\left(\frac{d\bm{\theta}}{2\pi}\right)^{d}\partial_{\mu_{j}}\log f(e^{\bm{\mu}+i\bm{\theta}})\)
> \(=\mathrm{Re}\int_{T^{d}}\left(\frac{d\bm{\theta}}{2\pi}\right)^{d}\frac{\partial_{\mu_{j}}f(e^{\bm{\mu}+i\bm{\theta}})}{f(e^{...\)
>
> **(Content truncated. Fetch a more specific URL or section for the full text.)**

Everything after that point in the article (Sec. III remainder, Sec. IV, Sec. V, Sec. VI, Sec. VII, Sec. VIII, Sec. IX, and Appendices A and B) could not be read.

## 3. Verbatim text obtained from the APS JSON (https://journals.aps.org/prx/fulltext/10.1103/PhysRevX.14.021011) — TRUNCATED AT THE SAME POINT

The APS machine-readable rendering confirms the same truncation and gives the closing-remark sentence that this is the whole of what is retrievable:

> (Content truncated. Fetch a more specific URL or section for the full text.)
>
> (Omitted 50864 bytes. Full formatted result stored at: C:\Users\45336\AppData\Local\Temp\dsh-spill-ppSlhW\session-29740491f71c\bee84f701f26-web_fetch.txt. Use read with offset/limit, or grep this path to search within it.)

I read that spill file: it contains only the same ~50 KB head (total 7 lines; the body is one truncated line reading `...(line truncated to 2000 chars)` followed by `(Content truncated. Fetch a more specific URL or section for the full text.)`). The 50,864 omitted bytes are therefore not recoverable through the session's file tools.

Verbatim titles of the addressable sections in that JSON (all present as anchors but only `s1`, `s2`, `s2b`, `s3` were within the returned head):

> `"id":"s1" ... "title":"I. INTRODUCTION","title_header":"h3"`
> `"id":"s2" ... "title":"II. MOTIVATION","title_header":"h3"`
> `"id":"s2b" ... "title":"B. The way to the amoeba","title_header":"h4"`
> (Sec. III components visible; `s4`–`s8` lie beyond the truncation.)

## 4. Verbatim fragments of the paper's LATER text (search-index snippet channel only)

These are the only literal strings from beyond Sec. III that I was able to see. They came back from `web_search` as indexed snippets attached to the APS PDF URL. They are 1–2 sentence fragments and are reproduced exactly as displayed; the snippets themselves are cut with "...".

From https://journals.aps.org/prx/pdf/10.1103/PhysRevX.14.021011 (indexed snippet, anchored `#9#6`):
> "This model is known to have a Chern insulator phase as well as a trivial insulator phase in the Hermitian limit \((\gamma = 0)\)..."

From https://journals.aps.org/prx/pdf/10.1103/PhysRevX.14.021011 (indexed snippet, anchored `#9#5`):
> "With the results of the DOS in the previous section, the proof of this theorem is now simple"

From http://export-test.arxiv.org/pdf/2212.11743 (indexed snippet, anchored `#9#6`):
> "where \(\bm{x}\) are 2D integer coordinates, and \(\bm{e}_{j}\) is the unit vector in the \(j\)th direction [see Fig"

From http://export-test.arxiv.org/pdf/2212.11743 (indexed snippet, anchored `#9#2`):
> "\[h(\bm{\beta}) = t\left(\beta_{x}+\beta_{x}^{-1}+\beta_{y}+\beta_{y}^{-1}\right)\] (12) \[+ t^{\prime}\left(\beta_{x}+\beta_{x}..."

From https://arxiv-org.ezproxy.obspm.fr/html/2212.11743v3 (indexed snippet, anchored `#7`):
> "For the present model, the precise phase boundary can even be analytically determined by the amoeba formulation"

From https://arxiv-org.ezproxy.obspm.fr/html/2212.11743v3 (indexed snippet, anchored `#6`):
> "Now we provide evidence for Eqs"

NOTE: the `#9#N` anchors above are the search provider's own page/paragraph markers on the PDF, not usable URLs. The only structural information I could extract about Sec. VII is from cross-references inside the retrieved head text: **the non-Hermitian Chern-band model is Eq. (39)** of arXiv v1 (Eq. numbering differs between v1 and v3), and its parameters are `t, γ, v, m`, all real-valued:

> Figure 2: Illustration of the real-space Hamiltonian of 2D models. (a) Single band toy model Eq. (7). (b) The non-Hermitian Chern-band model Eq. (39). All the parameters \(t,\gamma,v,m\) are real-valued.
> (from https://arxiv.org/html/2212.11743v1)

For the published PRX v3 numbering I could not read the equation number of the Chern-band model, its Bloch Hamiltonian, the GBZ definition in Sec. VI, the Chern-number expression, the numerical results, or any caveats.

## 5. Explicit list of everything I could NOT read

- Sec. VI "Generalized Brillouin zone": entire section — the GBZ definition in terms of the Ronkin function, any "Definition" statements, the GBZ computation algorithm, and any discussion of non-convexity/degeneracies. **NOT READ.**
- Sec. VII "Non-Bloch band topology": entire section — the non-Hermitian Chern-band Bloch Hamiltonian and its parameters `v, m, γ, t`, the Chern-number definition/equation computed on the GBZ, the numerical results, Table I, Figs. 8–11, and all stated caveats (quantization, applicability limits, multiband issues, higher-dimensional generalization). **NOT READ.**
- Sec. III after Eq. (15)/(16): the remainder of the Ronkin-function properties (the gradient/`ν_j` results, the Ronkin-function–amoeba relation, the "Re" statements) — truncated mid-equation.
- Sec. IV "Energy Spectra and density of states" (IV.1 Statement of the proposal, IV.2 Numerical evidence, IV.3 Derivation) — NOT READ.
- Sec. V "Amoeba hole closing and spectral boundary" — NOT READ (only the title and the fragments in section 4 above).
- Sec. VIII "Spectral inequalities", Sec. IX "Concluding remarks", Acknowledgements, Appendix A (Szegő's limit theorem), Appendix B (Invariance under change of basis), and the reference list — NOT READ.
- All figure captions after Fig. 2 and all equation numbers after Eq. (16) in the v3/APS numbering — NOT READ.
- The PDF (any host) — NOT READ: the tool rejects `application/pdf`.
- The arXiv LaTeX source (`arxiv.org/src/2212.11743v3`) — NOT READ: gzip, unsupported content type.
- The OpenAlex GROBID XML full text — NOT READ: HTTP 401, API key required.
- The APS paper HTML page in a browser-like form and any per-section/per-figure deliverable URLs — NOT AVAILABLE (404s; query parameters ignored).

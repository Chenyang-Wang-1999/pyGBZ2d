# arXiv:2609.23523v1 — full-text retrieval attempt across endpoints

Paper: "Universal Generalized Brillouin Zone Theory I: Review of the Spectral Approach"
Zeqi Xu, Jiangping Hu, Zhesen Yang — cond-mat.mes-hall, submitted 20 Sep 2026, 16 pages, 10 figures.
Verified to exist: `https://arxiv.org/abs/2609.23523v1` (HTTP 200).

## Verdict (short)

**No endpoint delivered text beyond Section II.** Sections III ("Boundary Matrix Equation", III.1, III.2)
and IV ("Exact Calculations for Finite Systems", IV.1, IV.2) appeared **only as table-of-contents
headings/anchor titles**, never as body text. No verbatim Section III or Section IV text exists in
anything this session received, so none is reported. This matches the independent capture already in
`notes_2609.23523_S1_S2.md` (lines 322–323 mark III and IV as not obtained).

## Exact truncation boundary (arXiv HTML rendering)

Every successful arXiv-HTML retrieval stopped at the *same* character position: the display equation
that follows the sentence

> "The hopping matrices $\{T_i\}$ also determine the real-space Hamiltonian under OBCs. For a chain of
> $L$ unit cells, with $L\geq p+q$ so that the boundary regions do not overlap, the eigenvalue equation reads"

i.e. at the unwinding of the first large matrix equation (the $nL\times nL$ block Toeplitz equation,
Eq. (13) in HTML order), inside **Section II**, after Eqs. (10), (11), (12) and the rank-deficiency
remark pointing to Sec. IV.2. The stall is caused by the enormous inline `\lxSVG@...` math fallback
markup of that matrix, which exhausts the fetch tool's character budget.

Last ~300 characters actually received (identical for every arXiv-HTML URL below):

```
\end{array}\right\]\left\[\begin{array}\[\]{c}\psi(1)\\\\[6.45831pt\] \psi(2)\\ \vdots\\ \hline\cr\\\[-6.45831pt\] \psi(p{+}1)\\\\[2.15277pt\] \vdots\\\[4.30554pt\] \psi(L{-}q)\\\\[2.15277pt\] \hline\cr\\\[-10.76385pt\] \vdots\\\[5.16663pt\] \psi(L{-}1)\\ \\\[-5.16663pt\] \psi(L)\\ \end{array}\right\]=E\left\[\begin{array}\[\]{c}\hbox to0pt{\vbox to0pt{\pgfpicture\makeatletter\hbox{\hskip 0.0pt\lower 0.0pt\hbox to0.0pt{\lxSVG@begingroup@{\_sc
```

## Per-endpoint results

| URL | Status | Depth reached | Tail (last ~300 chars) |
| --- | --- | --- | --- |
| https://arxiv.org/html/2609.23523v1 | HTTP 200, truncated | TOC (all 10 sections) + title/authors/abstract + I, I.1, I.2, II through Eq. (12) and the matrix of Eq. (13) | identical to block above |
| https://export.arxiv.org/html/2609.23523v1 | HTTP 200, truncated | byte-identical to previous, same stall point | identical to block above |
| https://arxiv.org/html/2609.23523 | HTTP 200, truncated | byte-identical, same stall point | identical to block above |
| https://arxiv.org/html/2609.23523v1?output=1 | HTTP 200, truncated | identical to base HTML (LaTeXML `output=1` changes nothing here) | identical to block above |
| https://arxiv.org/html/2609.23523v1#S3 and #S4 | HTTP 200, truncated | fragments are stripped server-side; identical full page, same stall point | identical to block above |
| https://ar5iv.labs.arxiv.org/html/2609.23523 | **Error**: cross-origin redirect to https://arxiv.org not followed automatically | retried against redirect target → https://arxiv.org/abs/2609.23523 (HTTP 200, abstract page only, no full text) | arXiv abs-page footer text |
| https://ar5iv.org/abs/2609.23523 | **Error**: cross-origin redirect to https://ar5iv.labs.arxiv.org not followed | retried https://ar5iv.labs.arxiv.org/abs/2609.23523 → **Error**: cross-origin redirect to https://arxiv.org. ar5iv has no HTML (and no abs page) for this paper | — |
| https://www.alphaxiv.org/pdf/2609.23523 | HTTP 200 | client-rendered shell only; page text ends at the title plus nav | `/ -  # Universal Generalized Brillouin Zone Theory I: Review of the Spectral Approach` |
| https://www.alphaxiv.org/abs/2609.23523 | HTTP 200 | metadata + abstract + "AI Overview: No overview yet"; no body text | footer/nav + similar-paper list; not paper text |
| https://arxiv.org/pdf/2609.23523v1 | **Error**: unsupported content type "application/pdf" (tool cannot ingest PDFs) | none | — |
| https://arxiv.org/e-print/2609.23523v1 (TeX source) | **Error**: unsupported content type "application/gzip" | none | — |
| https://www.emergentmind.com/papers/2609.23523 | HTTP 200 | abstract only | `[Abstract](https://arxiv.org/abs/2609.23523)[](https://arxiv.org/pdf/2609.23523)` |
| https://papers.cool/arxiv/2609.23523 | HTTP 200 | abstract + author/subject metadata only | `Designed by [kexue.fm](https://kexue.fm/) \| Powered by [kimi.ai](https://kimi.moonshot.cn/?ref=papers.cool)` |
| https://hjfy.top/arxiv/2609.23523 (and `?_x_output_type_...=embedded_pdf`) | HTTP 200 | client-rendered translation app shell | `     幻觉翻译 - arXiv 及文档翻译` |
| https://www.themoonlight.io/en/review/universal-generalized-brillouin-zone-theory-i-review-of-the-spectral-approach | HTTP 429 | Vercel bot checkpoint | `Vercel Security Checkpoint \| sin1::1790244930-SHpi3vRJhedhKC2oCzXWjEiRTXJAyx6q` |
| https://api.semanticscholar.org/graph/v1/paper/arXiv:2609.23523 | HTTP 429 | rate-limited | `Too Many Requests. Please wait and try again...` |
| https://arxiv-org.ezproxy.obspm.fr/html/2609.23523v1 | fetch failed (TypeError) | none; unreachable from this session | — |
| https://r.jina.ai/https://arxiv.org/html/2609.23523v1 and .../pdf/2609.23523v1 | fetch failed (TypeError) | none | — |
| https://www.textise.net/showText.aspx?strURL=... | fetch failed (TypeError) | none | — |
| https://www.arxiv-vanity.com/papers/2609.23523/ | fetch failed (TypeError) | none | — |
| https://synthical.com/article/2609.23523 | fetch failed (TypeError) | none | — |
| https://huggingface.co/papers/2609.23523 | fetch failed (TypeError) | none | — |

## Sections III / IV — appearance record

| Section | In body text? | Elsewhere |
| --- | --- | --- |
| III Boundary Matrix Equation | **No** | TOC entry + anchors `#S3`, `#S3.SS1`, `#S3.SS2`; cross-reference in I.2 |
| III.1 Boundary Matrix for a Finite Chain | **No** | TOC entry only |
| III.2 Questions in the Thermodynamic Limit | **No** | TOC entry only |
| IV Exact Calculations for Finite Systems | **No** | TOC entry + anchors `#S4`, `#S4.SS1`, `#S4.SS2`; cross-references in I.2 and II |
| IV.1 Example: Single-Band Model with Long-Range Hopping | **No** | TOC entry only |
| IV.2 Example: Non-Hermitian SSH Model with Rank Deficiency | **No** | TOC entry only, plus the Section II pointer to it |

## Why no endpoint can go deeper (diagnosis)

* The arXiv HTML is a single monolithic LaTeXML page; there is no section-split endpoint, and URL
  fragments (`#S3`, `#S4`) never reach the server.
* `?output=1` does not change the rendered document.
* The fetch tool's per-response character budget is consumed by the ~6 kB nav/TOC block plus the
  multi-kilobyte inline-SVG math fallback of Eq. (13); Section III therefore always falls past the cap.
* All PDF paths are barred at the transport level (`application/pdf` unsupported), and the TeX source
  is gzip (`application/gzip` unsupported), so neither can be converted to text in-session.
* Text-extraction proxies that would compress the markup (r.jina.ai, textise.net) are unreachable
  (TypeError: fetch failed), as is arxiv-vanity; ar5iv no longer covers recent arXiv IDs; alphaXiv,
  emergentmind, papers.cool and hjfy.top carry no body text.

## What would be needed to obtain Sections III–IV

1. A local PDF/TeX download outside the web-fetch path (curl/wget + pdftotext or `tar xzf`), or
2. A reachable HTML/text mirror that either splits the document by section or strips the inline-SVG
   math markup, or
3. A fetch tool that supports byte/character offsets or per-section extraction.

No content from Sections III or IV has been invented, inferred, or paraphrased here.

---

## Addendum (independent session, same tools): retrieval routes for arXiv:2212.11743

The same per-response character cap blocks the analogous read of Wang–Song–Wang, *Amoeba Formulation of
Non-Bloch Band Theory in Arbitrary Dimensions* (arXiv:2212.11743 → PRX **14**, 021011 (2024)). Confirmed
dead ends from a different session, recorded so nobody re-tries them:

| Route | Result |
|---|---|
| `https://arxiv.org/html/2212.11743v3`, `...?section=VI`, `.../?section=S6`, `.../`, `#S3`/`#S6`/`#S9` | HTTP 200 but **always the same window from the document start**, stalling mid-Sec. III (Ronkin definition) |
| `https://ar5iv.labs.arxiv.org/html/2212.11743` and `.../2212.11743v3` | same stall point |
| `https://journals.aps.org/prx/fulltext/10.1103/PhysRevX.14.021011` (+ `#s4`/`#s6`/`#s9`) | HTTP 200, returns the **fulltext JSON in document order**; anchors ignored; same stall point (mid-Sec. III) |
| `https://journals.aps.org/prx/abstract/...`, `https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevX.14.021011/fulltext` | `unsupported content type application/pdf` |
| `https://arxiv.org/pdf/2212.11743v3`, `https://arxiv.org/src/2212.11743` | `unsupported content type application/pdf` / `application/gzip` |
| `https://r.jina.ai/...`, `https://www.researchgate.net/...`, `https://www.semanticscholar.org/arxiv/2212.11743`, `https://phiscout.com/...` | fetch failed / 403 / 404 / DNS |
| `https://arxiv-org.ezproxy.obspm.fr/html/2212.11743v3` | cross-origin redirect to `login.ezproxy.obspm.fr` |
| `https://www.nature.com/articles/s42005-026-02546-2` (+ `.pdf`) | cross-origin redirect to `idp.nature.com` |
| `https://www.semanticscholar.org`, `https://www.alphaxiv.org/abs/2212.11743`, `https://www.emergentmind.com/papers/2212.11743` | no body text (alphaXiv overview is itself truncated inside Sec. IV) |

**Workaround that does yield primary-adjacent text:** papers by the same group / same line of work that
restate the 2212.11743 amoeba formulation, read via arXiv HTML (their Secs. I–II fit under the cap):

* arXiv:**2608.28577** (Gu, Fu, Hu, Wang) — https://arxiv.org/html/2608.28577v1 — restates the 2D amoeba
  formulation: `R_E(μx,μy) = ∫_{T²} dθx dθy/(2π)² log|f(e^{μx+iθx}, e^{μy+iθy})|`, `φ(E) = min_{μx,μy} R_E`,
  `ρ(E) = (1/2π)∇²φ(E)`; convexity of `R_E`; minimum = single point or finite plateau; plateau ⇔ E not in
  the OBC spectrum; single point ⇒ `|ψ_E(x,y)| ~ e^{μ̃x x + μ̃y y}`; 2D benchmark = DOS from Ronkin function
  vs ED on a square (L=130 with boundary disorder; and L=256, g=1e-4).
* arXiv:**2511.11349** (Kaneshiro & Peters, PRR 8, 013292 (2026)) — https://arxiv.org/html/2511.11349v2 —
  1D amoeba formulation with the **explicit Ronkin closed form**
  `R_σ(μ) = log|C_E| − p μ + Σ_{j=1}^{p+q} max(μ, μ_j)`, the fact that the minimum lies in `[μ_p, μ_{p+1}]`
  (= the 1D GBZ condition), and the statement that **implementation is limited to single-band systems**.
* arXiv:**2407.01296** (Xiong, Xing, Hu) — Sec. V "Amoeba formulation"; characterises 2212.11743's GBZ as
  "a d-dimensional object embedded in a 2d-dimensional space … without specifying the geometric information".
* Search snippets (unverified, page-tagged) from the PRX PDF mention an "explicit formula of the Ronkin
  function" applied to `det[E − h(β)] = a_{−M}(E)β^{−M} + …` and an extended formula
  `R_g(log|R|) = log|g(0)| + log|z₂/z₁| + 2log|z₃/z₂| + 3log|z₄/z₃| + …`.

None of the above is a substitute for Secs. IV–IX of 2212.11743 itself; it is recorded as secondary.

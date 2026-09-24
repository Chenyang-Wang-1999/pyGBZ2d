# Retrieval summary: Wang–Song–Wang, "Amoeba Formulation of Non-Bloch Band Theory in Arbitrary Dimensions"
arXiv:2212.11743v3 / Phys. Rev. X **14**, 021011 (2024). Open access CC-BY. 21 pages, 11 figures, 1 table, 62 refs.

Written by the retrieval agent. Companion files in the same directory:
`sec-IV-V.md`, `sec-VI-VII.md`, `sec-VIII-IX.md`, `sec-VIII-IX-FULL.md`,
`appendices-and-captions.md`, `sec-VII-VIII-appendices.md`,
`late-sections-via-search-index.md`, `derivative-papers-on-wsw.md`,
`late-equations-and-GBZ.md`, `ack-captions-appendices.md`.

---

## 1. Root cause of the whole problem

`web_fetch` returns only the **first ~50,000 characters** of any response body and has **no byte-range/offset support**. URL fragments (`#S7`) are stripped before the request, and the servers here **ignore query strings** (`?section=s7`, `?page=2`, `?format=json` on the APS fulltext all return byte-identical head output). The article's HTML is ~663,806 bytes (`content-length` reported by microlink), so every direct fetch stops inside Section II/III:

- `arxiv.org/html/2212.11743v3` → stops at Eq. (13) (end of Sec. II.2)
- `ar5iv.labs.arxiv.org/html/2212.11743` → stops in Sec. III at the amoeba definition, Eq. (14)
- `journals.aps.org/prx/fulltext/...` (machine-readable JSON, section anchors `s1`..`s8`) → stops inside Eq. (15)/(16)
- The harness's own "spill" dumps are **also** head-capped AND `read`/`grep` truncate a single physical line to 2000 chars, so omitted bytes are unrecoverable from disk.

## 2. Routes that WORKED (ranked)

1. **`api.microlink.io` server-side CSS-selector extraction** — THE bypass. It fetches and extracts on its server, so the 50 KB display cap never applies:
   `https://api.microlink.io/?url=https%3A%2F%2Farxiv.org%2Fhtml%2F2212.11743v3&meta=false&data.x.selector=%23S9&data.x.type=text`
   LaTeXML anchors: sections `#S1`..`#S9`, appendices `#A1`/`#A2`, figures `#S7.F9`…, acknowledgements `#S9.acknowledgements1`.
   **Caveat: anonymous daily quota (HTTP 429 `ERATE`) is now exhausted for this IP.** It intermittently admits requests and caches successes. Re-running `%23S7`, `%23S8`, `%23A1`, `%23A2`, `%23S5.F5`, `%23S7.F9`, `%23S7.F10`, `%23S7.F11` from a fresh IP should finish the paper.
2. **`journals.aps.org/prx/abstract/<DOI>`** — renders the article's **footnotes/endnotes verbatim** (notes 34, 38, 46, 47, 61) plus all 62 reference texts, the Popular Summary, and equation anchors (`#d44`, `#da11`, `#da12`).
3. **`web_search` phrase mining** — the index holds the publisher PDF (anchors like `#9#7` = page 9 paragraph 7) and the arXiv HTML (`#4`..`#8`) and returns **verbatim in-article passages** as the titles of its "Sources". It exposes only a fixed chunk neighbourhood (mostly PDF p.9), so it cannot be scrolled.
4. **arXiv/ar5iv HTML of *derivative* papers** — short enough to read past the truncation point; they restate WSW's formalism verbatim and cite its equation numbers.
5. `api.openalex.org/works/doi:10.1103/PhysRevX.14.021011` — metadata; confirms `has_content {pdf, grobid_xml}` with `content_urls` (see §4).
6. `journals.aps.org/prx/export/<DOI>?type=bibtex`, `harvest.aps.org/v2/journals/articles/<DOI>` — citation + funding metadata.
7. `www.w3.org/services/html2txt?url=...&noinlinerefs=on&nonums=on&start=<byte offset>` — reportedly supports a byte offset and worked for another paper today; here it returned HTTP 429 (Cloudflare) throughout.

## 3. Routes that FAILED (all verified)

`r.jina.ai` (fetch failed) · `harvest.aps.org/.../fulltext` (PDF → unsupported content type) · `content.openalex.org/works/W4394855923.grobid-xml` (**HTTP 401, API key required**) · `arxiv.org/src/2212.11743v3` (gzip → unsupported) · all PDF URLs (`application/pdf` unsupported) · `arxiv.org/html/2212.11743v3/{s7,s8,main,index}.html` (404 — no per-section files) · `journals.aps.org/prx/fulltext/<DOI>/{body,s6,s7}` (404) and `/supplemental/<DOI>` (404 — no supplemental) · `api.codetabs.com/v1/proxy` (522) · `api.allorigins.win/raw` (520/timeout) · `corsproxy.io` (403 keyless) · `cors.lol`/`cors.sh` (530) · `textise.net` (403 Cloudflare) · `eu-central-1.harvest...`/`arxiv-org.ezproxy.obspm.fr` (login redirect) · `export-test.arxiv.org`, `browse.arxiv.org`, `xxx.itp.ac.cn`, `www.arxiv-three.com` (DNS ENOTFOUND) · `core.ac.uk`, `colab.ws`, `europepmc.org`, `x-mol.com`, `baidu xueshu` (Cloudflare/CAPTCHA/403) · `semanticscholar.org/arxiv/2212.11743` (404) and its snippet API (429) · public SearXNG instances, `mojeek`, `brave`, `startpage`, `google` (403/fetch failed) · `docs.google.com/viewer` (fetch failed) · `arxiv-vanity`, `emergentmind`, `alphaxiv/pdf`, `huggingface.co/papers` (empty/JS shell/truncated) · `ar5iv.org/abs/...` and `dx.doi.org`/`link.aps.org` (cross-origin redirects not followed) · `subagent`/`subagent_fork` at this depth (refused: depth exceeds maxDepth).

## 4. The single best remaining lead

OpenAlex serves a **GROBID XML** full text (`<div><head>…</head><p>…</p></div>`, section-structured) at
`https://content.openalex.org/works/W4394855923.grobid-xml` — blocked here only by **HTTP 401 (API key required)**.
With a free OpenAlex key this should yield the *entire* article in one structured document. Second best: a shell with `curl -H "Range: bytes=…"` against `https://arxiv.org/html/2212.11743v3`.

## 5. Coverage ledger — what was recovered and what was not

### RECOVERED (verbatim)
- **Acknowledgements, complete**: "Acknowledgements. This work is supported by NSFC under Grant No. 12125405." — **no code-availability and no data-availability statement** in the arXiv v3 rendering; no APS metadata field either.
- **Section IX "Concluding remarks", complete** (4 paragraphs) — see `ack-captions-appendices.md`.
- **Section V** (complete per that agent) and **Section VI** (complete per that agent): GBZ definition (Eq. 41), Eq. (37) `rho(E)=0, E not in Lambda`, Eqs. (29)/(34) (original vs generalized Szegő), footnote 4 on exceptional points, Eq. (44) contour integral of the Green's function.
- **Footnotes 34, 38, 46, 47, 61**, verbatim (anti-aliasing of bulk randomness; zero-coupling limit; omitted boundary terms; nilpotent term at EPs; `sigma(e^{i theta}) = eps + e^{i theta}` heuristic for Eqs. (a11)/(a12)).
- **Figure captions 1, 2, 3, 4, 6, 7, 8** with geometries/system sizes: disk diameter L=140; square side L=130 with on-site random boundary potential uniform in [-0.5,0.5]; Fig. 7 fitting window L in [80,240]; Fig. 8 square L=130 and disk with L=400, points within distance 20 of the boundary discarded.
- **Sec. IV.2 "amoebic spectrum" paragraph** with the term's definition; Sec. IV.3 fragments (Toeplitz quasi-product bound `||T[sigma1 sigma2] - T[sigma1]T[sigma2]||_1 = O(L^{d-1})`, and a second inequality `||T[sigma]^{-1} - T[sigma^{-1}]||_1 <= ...` truncated).
- **Derivative-paper restatements** (clearly labelled as such): convexity of the Ronkin function and the single-point-vs-plateau criterion, `phi(E) = min_mu R_E(mu)`, `rho(E) = (1/2 pi) Delta phi(E)`, and an explicit multi-band failure statement.

### NOT RECOVERED
- **Section VII** body: Chern-band Hamiltonian/parameters, the Chern-number equation on the GBZ, numerical results, **Table I**, and all stated caveats. (Only: "This model is known to have a Chern insulator phase as well as a trivial insulator phase in the Hermitian limit (gamma = 0)…" and "For the present model, the precise phase boundary can even be analytically determined by the amoeba formulation".)
- **Section VIII** as complete inequalities: no §VIII equation number; `\tilde{\Phi}` still truncated at `\tilde{\Phi}_{\math`; only the two Toeplitz norm bounds and prose fragments ("An alternative proof of the spectral inequality is based on Eq…", "Therefore, a condition for Eq…", "Taking a circle |E| = R that surrounds the whole PBC spectrum, the average potential on this circle is…").
- **Section IV.1 and IV.3** bodies (the numbered statement of the theorem/proposal).
- **Appendix A and Appendix B** bodies; equations **(a11)/(a12)** themselves (only footnote 61's description of them).
- **Figure 5, 9, 10, 11** captions; the content of **Table I**.
- **The published PRX typesetting**: all quotes are the arXiv v3 (30 Apr 2024) rendering or the publisher PDF as seen through the search index; wording and equation numbers may differ slightly from PRX.

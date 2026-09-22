// ============================================================
//  zhaji-derived template for science and mathematics notes
//  Supports standalone lessons and combined books with visual heading levels
// ============================================================

// ---------- Shared state for book and standalone lesson layouts ----------
#let __is_book = state("__is_book", false)

// ---------- Fonts ----------
// Body: New Computer Modern, with Songti SC and PingFang SC fallbacks
#let font-text = ("New Computer Modern", "Songti SC", "PingFang SC")
// Headings: sans-serif font fallbacks
#let font-head = ("Heiti SC", "PingFang SC", "New Computer Modern")
// Mathematical expressions
#let font-math = ("New Computer Modern Math", "New Computer Modern")

#let thm-name = (
  "zh": "Theorem",
  "en": "Theorem",
)

#let def-name = (
  "zh": "Definition",
  "en": "Definition",
)

// ---------- Page colors ----------
#let c-accent = rgb("#222222")
#let c-remark = rgb("#777777")
#let c-blue   = rgb("#3b5f82")
#let c-amber  = rgb("#96704a")
#let c-emph   = rgb("#b02a2a") // Dark red emphasis color

// ---------- Top-level page and body layout ----------
#let note(
  title: "",               // Lesson or document title
  subtitle: "Lecture Notes",     // Book cover subtitle; none hides it
  author: "",              // Optional author
  date: auto,              // Cover date: auto for current month, none to hide, or custom text
  mode: "lesson",          // "lesson" for standalone notes or "book" for a collection
  font-head: font-head,
  font-math: font-math,
  font-text: font-text,
  font-size: 10.8pt,
  lang: "en",
  region: "us",
  first-line-indent: 2em,
  leading: 0.86em,
  par-spacing: 1.2em,
) = {
  // ---------- Callout blocks support content arguments and legacy named arguments ----------
  let hint(..args) = {
    let pos = args.pos()
    let named = args.named()
    let style = named.at("style", default: "gray")
    let (title, body) = if pos.len() >= 2 {
      (pos.at(0), pos.at(1))
    } else if pos.len() == 1 {
      (named.at("title", default: none), pos.at(0))
    } else {
      (none, [])
    }
    let colors = (gray: rgb("#777777"), blue: rgb("#526b84"), amber: rgb("#96704a"))
    let border-color = colors.at(style, default: rgb("#777777"))

    block(
      breakable: true,
      width: 100%,
      inset: (x: 0.75em, y: 0.45em),
      stroke: (left: 1pt + border-color),
    )[
      #set par(first-line-indent: 0em, justify: true)
      #if title != none {
        text(font: font-head, weight: "regular", fill: c-accent)[#title]
        h(0.6em)
      }
      #body
    ]
  }

  // ---------- Theorem and definition callouts ----------
  let thm(..args) = {
    let pos = args.pos()
    let named = args.named()
    let (title, body) = if pos.len() >= 2 {
      (pos.at(0), pos.at(1))
    } else if pos.len() == 1 {
      (named.at("title", default: none), pos.at(0))
    } else {
      (none, [])
    }
    hint(
      title: if title != none [#thm-name.at(lang) · #title] else [#thm-name.at(lang)],
      style: "blue",
      body,
    )
  }

  let def(..args) = {
    let pos = args.pos()
    let named = args.named()
    let (title, body) = if pos.len() >= 2 {
      (pos.at(0), pos.at(1))
    } else if pos.len() == 1 {
      (named.at("title", default: none), pos.at(0))
    } else {
      (none, [])
    }
    hint(
      title: if title != none [#def-name.at(lang) · #title] else [#def-name.at(lang)],
      style: "blue",
      body,
    )
  }

  let make(body) = {
    // Global body font and paragraph settings
    set text(font: font-text, size: font-size, lang: lang, region: region)
    set par(justify: true, leading: leading, first-line-indent: first-line-indent)
    set math.equation(numbering: none)
    show math.equation: set text(font: font-math)
    show math.equation.where(block: false): it => it
    show figure.caption: set align(left) 

    // Heading levels use visual styling without extra numbering
    // Level 1: Chapter or topic, with an accent rule below
    show heading.where(level: 1): it => block(width: 100%, above: 2.4em, below: calc.max(1.2em, par-spacing))[
      #set text(font: font-head, size: 20pt, weight: "bold", fill: c-accent)
      #it.body
      #v(0.35em)
      #line(length: 100%, stroke: 0.75pt + c-blue)
    ]

    // Level 2: Section, with a 3.5pt blue-gray bar and space above
    show heading.where(level: 2): it => block(width: 100%, above: 2.0em, below: calc.max(0.85em, leading))[
      #grid(
        columns: (auto, 1fr),
        gutter: 0.55em,
        align: (left + horizon, left + horizon),
        rect(width: 3.5pt, height: 1.15em, fill: c-blue, radius: 1pt),
        text(font: font-head, size: 15pt, weight: "bold", fill: c-accent)[#it.body],
      )
    ]

    // Level 3: Model or topic, with a square marker and 12.5pt text
    show heading.where(level: 3): it => block(above: 1.4em, below: calc.max(0.6em, leading))[
      // Raise the smaller square by 2pt to align it with the heading
      #text(fill: c-blue, size: 8.5pt, baseline: -2pt)[■]
      #h(0.45em)
      #text(font: font-head, size: 12.5pt, weight: "bold", fill: c-accent)[#it.body]
    ]

    // Level 4: Analysis step, with a dash marker and 11pt text
    show heading.where(level: 4): it => block(above: 1.0em, below: calc.max(0.45em, leading))[
      // Raise the smaller dash by 1pt to align it with the heading
      #text(fill: c-remark, size: 9pt, baseline: -1pt)[–]
      #h(0.35em)
      #text(font: font-head, size: 11pt, weight: "bold", fill: rgb("#444444"))[#it.body]
    ]


    if mode == "book" {
      __is_book.update(true)

      let doc-meta-title = if title != "" {
        title
      } else {
        "Lecture Notes"
      }

      set document(
        title: doc-meta-title,
        author: if author != "" { author } else { () },
      )

      // 1. Cover page without a header or footer
      set page(
        paper: "a4",
        margin: (x: 2.55cm, top: 2.2cm, bottom: 2.25cm),
        header: none,
        footer: none,
      )

      let cover-title = if title != "" {
        title
      } else {
        {
          let h1 = query(heading.where(level: 1))
          if h1.len() > 0 { h1.first().body } else { "Lecture Notes" }
        }
      }

      align(center + horizon)[
        #v(-2cm)
        #text(font: font-head, size: 28pt, weight: "bold")[#cover-title]
        #if subtitle != none and subtitle != "" [
          #v(1.2em)
          #text(font: font-text, size: 13.5pt, fill: c-remark)[#subtitle]
        ]
        #if author != "" [
          #v(2.5em)
          #text(font: font-text, size: 12pt, fill: c-accent)[#author]
        ]
        #v(5.5cm)
        #if date == auto [
          #text(font: font-text, size: 10pt, fill: c-remark)[
            #datetime.today().display("[month repr:long] [year]")
          ]
        ] else if date != none and date != "" [
          #text(font: font-text, size: 10pt, fill: c-remark)[#date]
        ]
      ]
      pagebreak()

      // 2. Contents page includes the first two heading levels
      outline(title: "Contents", depth: 2, indent: 1.5em)
      pagebreak()

      // 3. Body pages show the chapter title above and centered page numbers below
      set page(
        paper: "a4",
        margin: (x: 2.55cm, top: 2.2cm, bottom: 2.25cm),
        header: {
          let p = counter(page).get().first()
          let on-page = query(heading).filter(h => counter(page).at(h.location()).first() == p)
          let before-page = query(selector(heading).before(here()))
          let cur = if on-page.len() > 0 { on-page.first() } else if before-page.len() > 0 { before-page.last() } else { none }
          let head-text = if cur != none {
            cur.body
          } else if title != "" {
            title
          } else {
            "Lecture Notes"
          }
          set text(font: font-head, size: 8.5pt, fill: c-remark)
          align(left)[#head-text]
          v(-0.6em)
          line(length: 100%, stroke: 0.4pt + luma(70%))
        },
        footer: {
          align(center)[
            #set text(font: font-text, size: 8.5pt, fill: c-remark)
            #counter(page).display("1")
          ]
        },
      )
      counter(page).update(1)

      body
    } else {
      // Standalone lesson layout
      context {
        if __is_book.get() {
          // Included lessons inherit the book layout without resetting the page
          body
        } else {
          let hs2 = query(heading.where(level: 2))
          let hs1 = query(heading.where(level: 1))
          let running-title = if title != "" {
            title
          } else if hs2.len() > 0 {
            hs2.first().body
          } else if hs1.len() > 0 {
            hs1.first().body
          } else {
            "Lecture Notes"
          }

          set document(
            title: running-title,
            author: if author != "" { author } else { () },
          )

          set page(
            paper: "a4",
            margin: (x: 2.55cm, top: 2.2cm, bottom: 2.25cm),
            header: {
              set text(font: font-head, size: 8.5pt, fill: c-remark)
              align(left)[#running-title]
              v(-0.6em)
              line(length: 100%, stroke: 0.4pt + luma(70%))
            },
            footer: context {
              align(center)[
                #set text(font: font-text, size: 8.5pt, fill: c-remark)
                #counter(page).display("1")
              ]
            },
          )

          body
        }
      }
  }
  }

  (make: make, thm: thm, def: def, hint: hint)
}


// ---------- Emphasis and mathematical shorthand ----------
#let emph(body) = text(font: font-head, weight: "bold", fill: c-emph)[#body]
#let key(body) = box(inset: (x: 0.18em, y: 0.05em), radius: 2pt, fill: luma(92%))[#body]
#let qed = align(right)[$square$]

#let dd = math.dif                       // Differential operator d
#let pm = math.plus.minus                // Plus-minus sign
#let mp = math.minus.plus                // Minus-plus sign
#let R = math.bold(math.upright("R"))    // Real numbers
#let N = math.bold(math.upright("N"))
#let C = math.bold(math.upright("C"))
#let e = math.upright("e")               // Base of the natural logarithm
#let i = math.upright("i")
#let abs(x) = $|#x|$
#let norm(x) = $norm(#x)$
#let inner(a, b) = $angle.l #a, #b angle.r$

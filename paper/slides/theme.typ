// Academic theme based on the user's Touying University example.
#import "@preview/touying:0.7.4": *
#import themes.university: university-theme

#let ink = rgb("04364a")
#let blue = rgb("176b87")
#let teal = rgb("448c95")
#let muted = rgb("626975")

#let academic-theme = university-theme.with(
  aspect-ratio: "16-9",
  progress-bar: true,
  header: utils.display-current-heading(level: 2),
  header-right: none,
  footer-columns: (27%, 1fr, 18%),
  footer-a: self => self.info.author,
  footer-b: self => utils.display-current-heading(level: 1),
  footer-c: self => context [#utils.slide-counter.display() / #utils.last-slide-number],
  config-info(
    title: [From physical GBZ formulations\ to polynomial roots],
    short-title: [Numerical 2D GBZs],
    subtitle: [A numerical bridge for amoebic and strip GBZs],
    author: sys.inputs.at("author", default: "pyGBZ2d research project"),
    institution: sys.inputs.at("institution", default: none),
    date: datetime(year: 2026, month: 9, day: 7),
    logo: none,
  ),
  config-page(width: 960pt, height: 540pt,
    margin: (top: 67pt, bottom: 35pt, x: 38pt)),
  config-common(handout: true, new-section-slide-fn: none),
  config-methods(init: (self: none, body) => {
    set text(font: ("Libertinus Serif", "Cambria", "Times New Roman"), size: 22pt)
    set par(leading: .55em)
    set heading(numbering: "1.1")
    set math.equation(numbering: none)
    body
  }),
)

#let source(body) = text(size: 10pt, fill: muted, body)
#let two(left, right, widths: (1fr, 1fr)) = grid(columns: widths, gutter: 30pt, left, right)
#let label(body, color: ink) = block(above: 0pt, below: 7pt,
  text(size: 22pt, weight: "bold", fill: color, body))
#let note(body) = text(size: 17pt, fill: muted, body)
#let eq(body) = block(width: 100%, inset: (y: 4pt), body)
#let definition-block(title, body) = block(
  width: 100%, breakable: false,
  fill: rgb("f2f7f8"), stroke: (left: 2pt + ink),
  inset: 12pt,
  [#label(title)#body],
)

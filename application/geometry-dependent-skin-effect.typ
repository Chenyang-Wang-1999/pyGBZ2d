#import "style.typ": *  // modified from @preview/zhaji:0.1.0

#let nt = note(
  title: "Chern number calculation",
  author: "Chenyang Wang",
  font-head: "Arial",
	font-text: "Times New Roman",
	first-line-indent: 0em,
	lang: "en"
)

#show: nt.make
#show raw.where(block: true): it => block(it, fill: rgb("edf6fd"), inset: 5pt, breakable: false)
#show raw.where(block: false): it => highlight(it, fill: rgb("e0e0e0"), top-edge: 1.2em, bottom-edge: -.5em, extent: 0.2em)

#let mathbf(x) = $bold(upright(#x))$
#let rme = $upright(e)$
#let rmi = $upright(i)$

= Geometry-dependent skin effect



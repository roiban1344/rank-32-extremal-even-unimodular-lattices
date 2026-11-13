
#let font-serif = ("New Computer Modern", "IPAexMincho")
#let font-sans = ("DejaVu Sans", "IPAexGothic")

#let SL = math.op("SL")
#let GL = math.op("GL")
#let PSL = math.op("PSL")
#let PGL = math.op("PGL")
#let Aut = math.op("Aut")
#let End = math.op("End")
#let Hom = math.op("Hom")
#let Gal = math.op("Gal")
#let tr = math.op("tr")
#let Tr = math.op("Tr")
#let Nm = math.op("N")
#let rank = math.op("rank")
#let disc = math.op("disc")
#let Im = math.op("Im")
#let Ker = math.op("Ker")
#let OO = $cal(O)$
#let nm(x) = $(#x, #x)$

#let doc(title: none, subtitle: none, meta: none, body) = {
  set document(title: title)
  set page(paper: "a4", margin: (x: 22mm, top: 22mm, bottom: 24mm), numbering: "1")
  set text(font: font-serif, size: 10.5pt, lang: "ja", cjk-latin-spacing: auto)
  set par(justify: true, leading: 0.85em, first-line-indent: (amount: 1em, all: false), spacing: 1.2em)
  set heading(numbering: "1.1")
  show heading: set text(font: font-sans, weight: "bold")
  show heading.where(level: 1): it => { v(0.8em); it; v(0.3em) }
  set math.equation(numbering: "(1)", number-align: right)
  show math.equation: set text(font: ("New Computer Modern Math", "IPAexMincho"))
  show link: set text(fill: rgb("#1f4e9c"))
  show raw: set text(font: ("DejaVu Sans Mono", "IPAexGothic"), size: 0.92em)
  set table(stroke: 0.5pt + luma(160), inset: 5pt)
  show table: set text(size: 0.92em)
  set list(indent: 0.6em)
  set enum(indent: 0.6em)

  if title != none {
    align(center, block(width: 100%, {
      text(font: font-sans, size: 16pt, weight: "bold", title)
      if subtitle != none { linebreak(); v(0.2em); text(size: 11pt, subtitle) }
    }))
    v(0.4em)
  }
  if meta != none {
    block(width: 100%, inset: 8pt, fill: luma(245), radius: 3pt, text(size: 0.88em, meta))
    v(0.6em)
  }
  body
}

#let _box(label, color, body) = block(
  width: 100%, inset: (x: 9pt, y: 7pt), radius: 3pt,
  fill: color.lighten(88%), stroke: (left: 2.5pt + color),
  { text(font: font-sans, weight: "bold", size: 0.92em, fill: color.darken(20%), label); h(0.6em); body },
)
#let policy(body) = _box([Policy.], rgb("#3a6ea5"), body)
#let aside(body) = _box([Aside.], rgb("#4a8a4a"), body)
#let trnote(body) = _box([Translator's note.], rgb("#a0662a"), body)
#let caution(body) = _box([Caution.], rgb("#b03a3a"), body)

#let orig(p) = text(size: 0.8em, fill: luma(110), [[original #p]])

#let _thm-counter = counter("eeul-thm")
#let _thm(kind, name: none, body) = {
  _thm-counter.step()
  block(width: 100%, above: 1em, below: 1em, {
    text(font: font-sans, weight: "bold", { kind; " "; context _thm-counter.display() })
    if name != none { [ (#name)] }
    text(font: font-sans, weight: "bold", [.])
    h(0.5em)
    body
  })
}
#let theorem(name: none, body) = _thm([Theorem], name: name, body)
#let lemma(name: none, body) = _thm([Lemma], name: name, body)
#let proposition(name: none, body) = _thm([Proposition], name: name, body)
#let corollary(name: none, body) = _thm([Corollary], name: name, body)
#let definition(name: none, body) = _thm([Definition], name: name, body)
#let remark(name: none, body) = _thm([Remark], name: name, body)
#let example(name: none, body) = _thm([Example], name: name, body)
#let proof(body) = block(width: 100%, {
  text(font: font-sans, [Proof.]); h(0.5em); body; h(1fr); $square$
})

#let claim(kind, detail: none) = {
  let c = (
    "known": rgb("#2f6f3e"), "secondary": rgb("#7a6a2a"), "verified": rgb("#1f4e9c"),
    "sketch": rgb("#7a4a9c"), "conjecture": rgb("#9c4a7a"), "recalled": rgb("#b03a3a"),
    "discussion": rgb("#5a5a5a"), "open": rgb("#5a5a5a"),
  ).at(kind, default: luma(90))
  text(size: 0.85em, fill: c, {
    "["; kind
    if detail != none { [: #detail] }
    "]"
  })
}

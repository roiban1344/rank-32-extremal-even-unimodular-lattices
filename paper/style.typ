#import "style-base.typ": *

#let _lang = state("p32-lang", "en")

#let paper(lang: "en", title: none, subtitle: none, authors: none, date: none, body) = {
  show: doc.with(title: title, subtitle: subtitle, meta: none)
  set text(lang: lang)
  show figure.where(kind: table): set block(breakable: true)
  set heading(numbering: none)
  set math.equation(numbering: none)
  show figure.where(kind: table): set figure.caption(position: top)
  show figure.caption: set text(size: 0.92em)
  show figure.caption: set align(left)
  set figure(gap: 0.6em)
  _lang.update(lang)
  if authors != none or date != none {
    align(center, {
      if authors != none { authors; linebreak() }
      if date != none { text(size: 0.92em, date) }
    })
    v(0.8em)
  }
  body
}

#let small(body) = text(size: 0.9em, body)

#let _env(color, body) = block(
  width: 100%, inset: (x: 9pt, y: 7pt), radius: 2pt, breakable: true,
  fill: color.lighten(92%), stroke: (left: 2pt + color), body,
)
#let thm(body) = _env(rgb("#1f4e9c"), body)
#let fact(body) = _env(rgb("#2f6f3e"), body)
#let rem(body) = _env(luma(120), body)
#let proofbox(body) = block(width: 100%, inset: (left: 9pt, y: 2pt), stroke: (left: 0.6pt + luma(170)), breakable: true, [_Proof._ #body])
#let qed = h(1fr) + $square$
#let policy-box(body) = _env(rgb("#3a6ea5"), text(size: 0.95em, body))
#let abstract-box(body) = block(width: 100%, inset: (x: 18pt, y: 8pt), breakable: true, {
  align(center, text(font: font-sans, weight: "bold", context if _lang.get() == "ja" { [Abstract] } else { [Abstract] }))
  v(0.2em)
  set par(first-line-indent: 0pt)
  body
})

#let tagged(t, eq) = block(width: 100%, grid(columns: (1fr, auto), align: (center + horizon, right + horizon), eq, [(#t)]))

#let ptable(columns: 1, align: auto, caption: none, ..cells) = figure(
  kind: table, numbering: none, caption: caption,
  table(columns: columns, align: align, ..cells),
)

#let refs-list(..items) = {
  set par(first-line-indent: 0pt, hanging-indent: 1.5em)
  for it in items.pos() { block(above: 0.5em, below: 0.5em, it) }
}

#let pending(body) = highlight(fill: rgb("#fff1a8"), body)

#let todo(body) = block(width: 100%, inset: (x: 8pt, y: 5pt), radius: 2pt,
  fill: rgb("#f3f3f3"), stroke: (left: 2pt + luma(150)),
  text(size: 0.88em, fill: luma(70), [#strong[To be written.] #body]))

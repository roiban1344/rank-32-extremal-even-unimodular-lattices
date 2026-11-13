#import "style.typ": *
#import "computations.typ": calcs
#import "computations-p2.typ": calcs2
#import "computations-ja.typ": calc-ja, calc2-ja, calc-time-ja, release-content-ja

#let proofbox(body) = block(width: 100%, inset: (left: 9pt, y: 2pt), stroke: (left: 0.6pt + luma(170)), breakable: true,
  [#text(font: font-sans, weight: "bold", size: 0.95em)[証明.] #body])

#let abstract-box(body) = block(width: 100%, inset: (x: 18pt, y: 8pt), breakable: true, {
  align(center, text(font: font-sans, weight: "bold")[要旨])
  v(0.2em)
  set par(first-line-indent: 0pt)
  body
})

#let ja-rules(body) = {
  show emph: set text(font: ("New Computer Modern", "IPAexGothic"))
  show strong: set text(font: ("New Computer Modern", "IPAexGothic"))
  show table.cell: set par(justify: false)
  body
}

#let calc-table-ja() = {
  block(above: 0.8em, below: 0.4em, text(size: 0.92em, [表 11. 計算とコード。本文の \[C1\] などの参照は、この表の行を指す。]))
  table(columns: (auto, 1fr, auto), align: (left, left, right),
    table.header([ID], [計算とコード], [CPU 時間]),
    ..calcs.map(c => (
      [#c.id #label("calc-" + c.id)],
      {
        calc-ja.at(c.id, default: c.en)
        linebreak()
        text(size: 0.8em, c.files.map(f => raw(f.replace("research record: ", "研究記録: "))).join(linebreak()))
      },
      calc-time-ja.at(c.id, default: c.time),
    )).flatten()
  )
}

#let calc2-table-ja(caption: none) = {
  block(above: 0.8em, below: 0.4em, text(size: 0.92em, caption))
  table(columns: (auto, 1fr, auto), align: (left, left, right),
    table.header([ID], [計算とコード], [CPU 時間]),
    ..calcs2.map(c => (
      [#c.id #label("calc2-" + c.id)],
      {
        calc2-ja.at(c.id, default: c.en)
        linebreak()
        text(size: 0.8em, c.files.map(f => raw(f.replace("research record: ", "研究記録: "))).join(linebreak()))
      },
      calc-time-ja.at(c.id, default: c.time),
    )).flatten()
  )
}

#let release-content(f, markup: false) = {
  if f.name in release-content-ja { release-content-ja.at(f.name) }
  else if markup { eval(f.content, mode: "markup") } else { f.content }
}

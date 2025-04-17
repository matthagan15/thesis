#import "macros.typ": *
#import "conf.typ": *

#let ket(psi) = $lr(|#psi angle.r)$
#let bra(psi) = $lr(angle.l #psi|)$
#let ketbra(a, b) = $|#a angle.r angle.l #b|$
#let braket(a, b) = $angle.l #a|#b angle.r$
#let tp = $times.circle$
#let id = $bb(1)$
#let dmd = $diamond.medium$
// Common objects
#let hilb = $cal(H)$
#let partfun = $cal(Z)$
#let identity = $bb(1)$
#let gue = $"GUE"$
#let sinc = math.op("sinc")
#let hermMathOp = math.op("Herm")
#let im = math.op("Im")
#let diag = math.op("diag")
#let herm(x) = $hermMathOp parens(#x)$


#import "@preview/ctheorems:1.1.3": *
// #let lemma_og = thmbox("lemma", "Lemma", stroke: 1pt, bodyfmt: x => text(x, style: "italic"))
#let proof = thmproof("proof", "Proof", inset: (x: 0cm))
// #let theorem_og = thmbox("theorem", "Theorem", stroke: 1pt, bodyfmt: x => text(x, style: "italic"))
// #let definition_og = thmbox("definition", "Definition", stroke: 1pt, bodyfmt: x => text(x, style: "italic"))

#let lemma = thmbox("lemma", "Lemma", stroke: 1pt, bodyfmt: x => text(x, style: "italic"), fill: rgb("e8887377"))
#let theorem = thmbox(
    "theorem",
    "Theorem",
    stroke: 1pt,
    bodyfmt: x => text(x, style: "italic"),
    fill: rgb("#c8f6ad"),
)
#let definition = thmbox(
    "definition",
    "Definition",
    stroke: 1pt,
    bodyfmt: x => text(x, style: "italic"),
    fill: rgb("#62b6cb44"),
)

#let todo = x => { text([TODO: #x], fill: red, weight: "bold") }
#set math.equation(number-align: bottom)

#let q0 = $1 / (1 + e^(-beta gamma))$
#let q1 = $e^(-beta gamma) / (1 + e^(-beta gamma))$

#show: thmrules.with(qed-symbol: $square$)

#counter(heading).update(0)
#set heading(numbering: "A", supplement: "Appendix")
#show heading.where(level: 1): it => {
    let number = if it.numbering != none {
        context counter(heading).display(it.numbering)
        // h(1em)
    }
    pagebreak()
    v(3.5cm)
    block(text("Appendix ", size: 22pt) + text(number, size: 22pt))
    v(1.0cm)
    block(text(it.body, size: 25pt))
    v(14mm)
}
#set page(
    numbering: "1",
    margin: (left: 35mm, top: 20mm, bottom: 20mm, right: 20mm),
    footer: context if is-chapter-page() {
        align(center)[
            #counter(page).display("1") #h(15mm)
        ]
    } else { [] },
    header: context if is-chapter-page() {
        []
    } else {
        align(right)[
            #smallcaps([Appendix ] + hydra(1)) #h(1fr) #counter(page).display("1")
        ]
    },
)
= Trotter Bounds

= Haar Integrals <sec_appendix_haar>

#lemma([Sinc Function Bounds])[
    For $sinc^2(x t/2)$ and $delta_min$ as defined in @eq_delta_min_def, we will make significant use of the following Bounds
    $
        |x| >= delta_"min" ==> sinc^2 ( frac(x t, 2) ) &<= frac(4, delta_"min"^2 t^2) #label("eq:sinc_upper_bound") \
  |x| <= frac(sqrt(2), t) ==> sinc^2(frac(x t, 2) ) &>= 1 - frac(|x|^2 t^2, 2).
    $
] <lem_sinc_poly_approx>

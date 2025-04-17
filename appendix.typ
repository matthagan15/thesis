#import "macros.typ": *
#import "conf.typ": *

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

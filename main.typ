#import "conf.typ": ut-thesis, is-chapter-page
// #import("conf.typ")
#import "@preview/hydra:0.6.1": hydra
#include "macros.typ"

// Example Usage
#show: ut-thesis.with(
    title: "Preparing Thermal States on a Digital Quantum Computer",
    author: "Matthew Hagan",
    degree: "Doctor of Philosophy",
    department: "Physics",
)


#include "intro.typ"
#include "composite.typ"
#include "thermal_state.typ"
#include "conclusion.typ"

#show heading.where(level: 1): it => {
    let number = if it.numbering != none {
        context counter(heading).display(it.numbering)
    }
    pagebreak()
    v(1cm)
    block(text(it.body, size: 25pt))
    v(1cm)
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
            #smallcaps(hydra(1)) #h(1fr) #counter(page).display("1")
        ]
    },
)

#include "appendix.typ"

#bibliography("references.bib", style: "american-physics-society")

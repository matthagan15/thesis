#import "conf.typ": ut-thesis
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
#include "appendix.typ"

#bibliography("references.bib")

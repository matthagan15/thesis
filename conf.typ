// University of Toronto Thesis Typst Template
#import "@preview/hydra:0.6.1": hydra

/// returns bool
#let is-chapter-page() = {
    // all chapter headings
    let chapters = query(heading.where(level: 1))
    // return whether one of the chapter headings is on the current page
    chapters.any(c => c.location().page() == here().page())
}

/// Main Document Structure
#let ut-thesis(
    title: "[Thesis Title]",
    author: "[Author Name]",
    degree: "Doctor of Philosophy",
    department: none,
    graduation-year: datetime.today().year(),
    body,
) = {
    // Page setup
    set page(
        paper: "us-letter",
        margin: (left: 32mm, top: 20mm, bottom: 20mm, right: 20mm),
    )
    set text(
        font: "New Computer Modern",
        top-edge: 0.7em,
        bottom-edge: -0.3em,
        size: 11pt,
    )

    // Title Page (no numbering, centered)
    page(
        margin: (x: 0cm, y: 0cm),
        footer: none,
        align(center)[
            #v(30mm)
            #text(size: 12pt)[#smallcaps[#title]]

            #v(4cm)
            by

            #v(4cm)
            #text(size: 12pt)[#author]

            #v(5cm)
            A thesis submitted in conformity with the requirements \
            for the degree of #degree\
            Department of #department \
            University of Toronto

            #v(3cm)
            © Copyright by #author #graduation-year \
        ],
    )

    set page(
        numbering: "i",
        number-align: center,
    )

    page[
        #align(center)[
            #v(1em)
            #text(size: 12pt)[
                #title\
                #v(0.075cm)
                #author\
                #degree\
                #v(0.075cm)
                Department of Physics \
                University of Toronto\
                #graduation-year
            ] \
            #text(size: 14pt, weight: "bold")[Abstract]
        ]
        #set par(leading: 1.5em)
        The simulation of quantum systems is the most promising candidate for economic advantages of large scale, fault tolerant quantum computers over classical computers. The development of Markov Chain Monte Carlo techniques is one of the central tools for studying quantum systems on classical computers. We adapt modern techniques from classical Monte Carlo algorithms, specifically the Hamiltonian Monte Carlo (HMC) algorithm, to an algorithm for preparing thermal states of quantum systems on digital quantum computers. The simulation of time dynamics is a crucial subroutine for preparing thermal states using HMC, and to this end we provide new techniques for combining existing product formulas to achieve lower costs for overall simulations. These methods provide an avenue for extending the Repeated Interactions framework for studying quantum thermodynamics to arbitrary systems.
    ]

    page[
        #place(horizon + right, [_This thesis is dedicated to \ my brother JT and my sisters \ Brittany and Veronica._])
    ]

    page[
        #quote(
            block: true,
            attribution: [William Faulkner, The Sound and the Fury],
        )[_"...I give you the mausoleum of all hope and desire...I give it to you not that you may remember time, but that you might forget it now and then for a moment and not spend all of your breath trying to conquer it. Because no battle is ever won he said. They are not even fought. The field only reveals to man his own folly and despair, and victory is an illusion of philosophers and fools."_]
    ]


    page[
        #set par(leading: 1.5em, justify: true, first-line-indent: 5mm)
        #align(
            center,
            text(size: 14pt, weight: "bold", [Acknowledgements]),
        )
        There are many people I would like to thank for helping me through my Ph.D. First and foremost is my family, including but not certainly not limited to, my parents, Frank, Rhonda, Brittany, Veronica, JT, and my Grandfather. You all serve as my compass, and I cannot remember how many times I have had a difficult problem or decision to make and I talk to all of y'all to figure out what I should do. You all have taught me what the right priorities to have in life are and, most importantly, how to stay grounded while pursuing my dreams.

        I would like to also thank all of the friends throughout my graduate school journey. Starting with the cohort at the University of Washington who made our first year hallway sheer entertainment, despite the constant pressure of assignments and teaching. A special thanks is in order to my friends in Seattle, JT, Bailey, Diego, Caitlin, and Ann, for all the talks, hikes, parties, and Super Smash Brothers. I would also like to thank the friends I have made in Toronto for helping me adjust to moving to a new city and new country. To Nick, my drumming and life mentor. To Aaron, Luke, Robyn, Julian, Griffin, Asenia, Deepanshu, Andrija, and Joscelyn for making the neighborhood a fun place to be. To Alisha, for making my last year in graduate school the best it could be. And to my abhi, Burak. I have routinely wondered how someone raised on the other side of the world could value the same things in life that I do, the same way I do. You have been a (literal) constant in my life and have picked me up from rock bottoms with nothing but a sense of brotherhood and compassion. I hope I have returned the favor even half as much.

        Before I even applied to graduate school I did not think I was capable of even being accepted for, and much less completing, a Ph.D, so a special thanks to Maya Sathaye and Ibrahim Cisse for believing in me strongly enough to change my mind. The ideas for this thesis began in the Spring of 2016, when my roommate Michael Traub introduced me to an algorithm known as Hamiltonian Monte Carlo. There are many other things to thank you for, but a special thanks for planting the seeds of an idea that became this thesis.

        There are a number of people at the University of Toronto I would like to thank. To the members of my research group, thank you for putting up with my hare-brained ideas and claims that every algorithm is ultimately a thermal state preparation routine. My utmost appreciation for all of you, Juan Castañeda, Matt Pocrnic, Serene Shum, Sophia Simon, Abhinav Muraleedharan, Hokiat Lim, and Raghav Banka. I would also like to thank Dvira Segal, for somehow finding the time meet and discuss issues that sometimes came up while I developed these ideas further and to Stephen Julian for providing valuable feedback during our meetings.

        And lastly I would like to thank Nathan. I am incredibly grateful I decided to skip the remaining "dog and pony" sessions to chat science with someone who had seriously presented on simulating quantum mechanics in Minecraft. Thank you for your commitment to teaching me your craft, for pushing me when I needed it, for giving me the freedom to pursue my own crazy ideas, and for supporting these ambitions however possible.
    ]

    page[
        #outline(indent: 7.5mm)
    ]

    page[
        #outline(title: "List of Figures", target: figure.where(kind: image))
    ]

    page[
        #outline(title: "List of Tables", target: figure.where(kind: table))
    ]


    counter(page).update(1)
    set page(
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
                #smallcaps([Chapter ] + hydra(1)) #h(1fr) #counter(page).display("1")
            ]
        },
    )
    set par(leading: 1em, justify: true, first-line-indent: 5mm)

    // Render the body of the document
    set heading(numbering: (..nums) => nums.pos().map(str).join("."))
    show heading.where(level: 1): it => {
        let number = if it.numbering != none {
            context counter(heading).display(it.numbering)
            // h(1em)
        }
        pagebreak()
        v(3.5cm)
        block(text("Chapter ", size: 22pt) + text(number, size: 22pt))
        v(1.0cm)
        block(text(it.body, size: 25pt))
        v(14mm)
    }
    show heading.where(level: 2): it => {
        let number = if it.numbering != none {
            context counter(heading).display(it.numbering)
        }
        v(0.25cm)
        block(text(number + h(0.2cm) + it.body, size: 16pt))
        v(0.25cm)
    }
    // set math.equation(numbering: "(1.1)")
    set math.equation(
        numbering: it => {
            let count = counter(heading.where(level: 1)).at(here()).first()
            if count > 0 {
                numbering("(1.1)", count, it)
            } else {
                numbering("(1)", it)
            }
        },
    )
    body
}

#let todo = x => { text([TODO: #x], fill: red, weight: "bold") }

#let ket(psi) = $lr(|#psi angle.r)$
#let bra(psi) = $lr(angle.l #psi|)$
#let ketbra(a, b) = $|#a angle.r angle.l #b|$
#let braket(a, b) = $angle.l #a|#b angle.r$
#let bracket(a, b, c) = $angle.l #a|#b|#c angle.r$
#let tp = $times.circle$
#let id = $bb(1)$
#let dmd = $diamond.medium$
#let hilb = $cal(H)$
#let partfun = $cal(Z)$
#let identity = $bb(1)$
#let gue = $"GUE"$
#let sinc = math.op("sinc")
#let hermMathOp = math.op("Herm")
#let im = math.op("Im")
#let diag = math.op("diag")
#let herm(x) = $hermMathOp parens(#x)$
#let tpose = sym.top
#let on = $"on"$


#import "@preview/ctheorems:1.1.3": *
#let proof = thmproof("proof", "Proof", inset: (x: 0cm))

#let lemma = thmbox(
    "lemma",
    "Lemma",
    stroke: 1pt,
    base_level: 1,
    bodyfmt: x => text(x, style: "italic"),
    // fill: rgb("e8887377"),
)

#let theorem = thmbox(
    "lemma",
    "Theorem",
    stroke: 1pt,
    base_level: 1,
    bodyfmt: x => text(x, style: "italic"),
    // fill: rgb("#c8f6ad"),
)

#let definition = thmbox(
    "lemma",
    "Definition",
    stroke: 1pt,
    base_level: 1,
    bodyfmt: x => text(x, style: "italic"),
    // fill: rgb("#62b6cb44"),
)

#let corollary = thmbox(
    "lemma",
    "Corollary",
    stroke: 1pt,
    bodyfmt: x => text(x, style: "italic"),
    // fill: rgb("#c4c67d"),
)

// University of Toronto Thesis Typst Template
#import "@preview/hydra:0.6.1": hydra
#set text(font: "New Computer Modern")

/// returns bool
#let is-chapter-page() = {
    // all chapter headings
    let chapters = query(heading.where(level: 1))
    // return whether one of the chapter headings is on the current page
    chapters.any(c => c.location().page() == here().page())
}

/// Main Document Structure
#let ut-thesis(
    title: "Thesis Title",
    author: "Author Name",
    degree: "Doctor of Philosophy",
    department: "Graduate Department",
    graduation-year: datetime.today().year(),
    body,
) = {
    // Page setup
    set page(
        paper: "us-letter",
        margin: (left: 32mm, top: 20mm, bottom: 20mm, right: 20mm),
    )
    set text(font: "New Computer Modern")

    // Title Page (no numbering, centered)
    page(
        // Approximate positioning based on guidelines

        align(center)[
            #v(30mm)
            #text(size: 12pt)[#smallcaps[#title]]

            #v(4cm)
            by

            #v(4cm)
            #text(size: 12pt)[#author]

            #v(5cm)
            A thesis submitted in conformity with the requirements \
            for the degree of Doctor of Philosophy\
            Department of Physics \
            University of Toronto

            #v(3cm)
            © Copyright by #author #graduation-year
        ],
    )

    // Preliminary pages with Roman numerals
    set page(
        numbering: "i",
        number-align: center,
    )

    // Abstract (double-spaced, max 350 words for doctoral thesis)
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

        #set par(leading: 2em) // Double-spaced
        // Placeholder for abstract content
        // Actual content should be added by the user
        #lorem(100)
    ]

    page[
        #place(horizon + right, "This thesis is dedicated to...")
    ]

    page[
        #align(
            center,
            text(size: 14pt, weight: "bold", [Acknowledgements]),
        )
    ]

    page[
        #outline()
    ]

    // Optional lists (user can customize/remove as needed)
    page[
        #outline(title: "List of Figures", target: figure.where(kind: image))
    ]

    page[
        #outline(title: "List of Tables", target: figure.where(kind: table))
    ]


    // Main document body with Arabic numerals
    counter(page).update(1)
    set page(
        numbering: "1",
        margin: (left: 35mm, top: 20mm, bottom: 20mm, right: 35mm),
        footer: context if is-chapter-page() {
            align(center)[
                #counter(page).display("1")
            ]
        } else { [] },
        header: context if is-chapter-page() {
            []
        } else {
            align(right)[
                CHAPTER #upper(hydra(1)) #h(1fr) #counter(page).display("1")
                // CHAPTER #hydra(1, display: (hydra_ctx, hydra_content) => hydra_content) #h(1fr) #counter(page).display("1")
            ]
        },
    )
    // set page(number-align: top + right)
    set text(
        top-edge: 0.7em,
        bottom-edge: -0.3em,
        font: "New Computer Modern",
        size: 11pt,
    )
    set par(leading: 1em, justify: true, first-line-indent: 7mm)

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
    body
}

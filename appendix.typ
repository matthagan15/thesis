#import "conf.typ": *

#import "@preview/ctheorems:1.1.3": *


#let q0 = $1 / (1 + e^(-beta gamma))$
#let q1 = $e^(-beta gamma) / (1 + e^(-beta gamma))$

#show: thmrules.with(qed-symbol: $square$)

#counter(heading).update(0)
#set heading(numbering: "A.1", supplement: "Appendix")
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

#set math.equation(number-align: bottom)

= Sinc Approximation <sec_appendix_sinc>

#lemma([Sinc Function Bounds])[
    For $sinc^2(x t / 2)$ and $delta_min$ as defined in @eq_delta_min_def, we will make significant use of the following Bounds
    $
        |x| >= delta_"min" ==> sinc^2 ( frac(x t, 2) ) &<= frac(4, delta_"min"^2 t^2) \
        |x| <= frac(sqrt(2), t) ==> sinc^2(frac(x t, 2) ) &>= 1 - frac(|x|^2 t^2, 2).
    $ <eq_sinc_bounds>
] <lem_sinc_poly_approx>
#proof()[
    The first inequality is rather trivial
    $
        "sinc"^2 ((x t) / 2) = (sin^2((x t) / 2) ) / ((x t) / 2)^2 <= 4 / (x^2 t^2) <= frac(4, delta_"min"^2 t^2).
    $

    The second involves a Taylor Series for $"sinc"^2$, which we compute using the expression of $"sinc"$ as $"sinc"((x t) / 2) = ("sin" (x t) / 2) / ((x t) / 2) = integral_0^1 "cos"(s (x t) / 2) d s$. The first two derivatives can then be computed easily
    $
        frac(d "sinc"^2((x t) / 2), d x) &= - t integral_0^1 "sin"(s x) s " "d s integral_0^1 "cos"(s x) d s \
        frac(d^2 "sinc"^2((x t) / 2), d x^2) &= -t^2 / 2 integral_0^1 "cos"(s x) s^2 d s integral_0^1 "cos"(s x) d s + t^2 / 2 integral_0^1 "sin"(s x) s d s integral_0^1 "sin"(s x) s d s.
    $

    We can evaluate each of these derivatives about the origin using continuity of the derivatives along with the limits $lim_(x -> 0) "cos"(s x) = 1$ and $lim_(x -> 0) "sin"(s x) = 0$. We can now compute the mean-value version Taylor series as

    $
        "sinc"^2 (frac(x t, 2)) &= "sinc"^2(0) + x frac(d, d x) "sinc"^2 (frac(x t, 2)) |_(x = 0) + frac(x^2, 2!) frac(d^2, d x^2) "sinc"^2 (frac(x t, 2)) |_(x = x_star),
    $

    where $x_star in [0,1]$.
    Plugging in $"sinc"^2(0) = 1$ and $frac(d"sinc"^2(x t / 2), d x)|_(x = 0) = 0$ then yields
    $
        abs("sinc"^2((x t) / 2) - 1) = frac(|x|^2, 2) abs((d^2"sinc"^2(x t / 2)) / (d x^2) |_(x = x_star)).
    $
    We make use of the rather simplistic bound

    $
        abs(frac(d^2"sinc"^2(s x t / 2), d x^2)|_(x = x_star)) &<= t^2 / 2 abs(integral_0^1 "cos"(s x_star t / 2) s^2 d s integral_0^1 "cos"(s x_star t / 2) d s) + t^2 / 2 abs(integral_0^1 "sin"(s x_star t / 2) s d s integral_0^1 "sin"(s x_star t / 2) s d s) \
        &<= t^2 / 2 integral_0^1 abs("cos"(s x_star t / 2)) s^2 d s integral_0^1 abs("cos"(s x_star t / 2 )) d s + t^2 / 2 (integral_0^1 abs("sin"(s x_star t / 2)) |s| d s)^2 \
        &<= t^2 / 2 integral_0^1 s^2 d s + t^2 / 2 (integral_0^1 s d s)^2 \
        &<= t^2.
    $

    This yields the final inequality $|"sinc"^2((x t) / 2 ) - 1| <= frac(|x|^2 t^2, 2)$ which yields @eq_sinc_bounds.
]

= Random Interaction Model <sec_appendix_random_interaction>
In this appendix we show that the random interaction model $G$ satisfies the conditions needed in @thm_tsp_first_order_phi and .
#lemma()[
    Let $G = U_G Lambda_G U_G^dagger$ be given as described in @sec_tsp_prelims, specifically let $lambda_G(i) = bra(i) Lambda_G ket(i)$ denote the $i^"th"$ eigenvalue. Then the following expectation values hold
    $
        EE_G [G] = 0 " and " EE_(Lambda_G) [lambda_G (i) lambda_G (j)] = delta_(i, j).
    $
] <lem_tsp_interaction_expectations>
#proof()[
    The random interaction $Lambda_G$ is a collection of Pauli $Z$ strings and an overall phase of $plus.minus 1$ as
    $
        Lambda_G = (-1)^(z_0) Z_1^(z_1) tp ... tp Z_n^(z_n),
    $
    where $"Pr"[z_i = 0] = "Pr"[z_i = 1] = 1 / 2$. We first show that $EE_(Lambda_G) [Lambda_G] = 0$, which ultimately comes from the expectation over the phase. This can be computed as the $z_i$ are independent and the expectation therefore factors
    $
        EE_(Lambda_G) [Lambda_G] &= EE_(z_0) (-1)^(z_0) product_(i=1)^n EE_(z_i) tp Z_i^(z_i) \
        &= 1 / 2 dot (+1) product_(i=1)^n EE_(z_i) tp Z_i^(z_i) + 1 / 2 dot (-1) product_(i=1)^n EE_(z_i) tp Z_i^(z_i) \
        &= 0.
    $
    This implies that the overall interaction has a zero first moment via $EE_G [G] = EE_(U_G) [ U EE_(Lambda_G )[Lambda_G] U_G^dagger]$.

    Next we need to show that the eigenvalues are independent and have variance 1. First we let $z dot k = z_0 + z_1 k_1 + ... + z_n k_n$. First we show how this gives the eigenvalues of $Lambda_G$ as
    $
        Lambda_G ket(i) &= (-1)^(z_0) Z_1^(z_1) tp ... tp Z_n^(z_n) ket(i_1) tp ... tp ket(i_n) \
        &= (-1)^(z_0) (-1)^(z_1) ... (-1)^(z_n) ket(i_1) tp ... tp ket(i_n) \
        &= (-1)^(z dot i) ket(i).
    $
    This allows us to compute the covariance as
    $
        EE_(Lambda_G) bra(i) Lambda_G ket(i) bra(j) Lambda_G ket(j) &= EE_(Lambda_G) (-1)^(z dot i) (-1)^(z dot j) \
        &= EE_(Lambda_G) (-1)^(z dot (i + j)) \
        &= EE_(z_0) (-1)^(2 z_0) product_(k = 1)^n EE_(z_k) (-1)^(z_k (i_k + j_k)) \
        &= product_(k = 1)^n EE_(z_k) (-1)^(z_k (i_k + j_k)).
    $
    Now we just need to compute a single one:
    $
        EE_(z_k) (-1)^(z_k (i_k + j_k)) &= 1 / 2 dot (1) + 1 / 2 dot (-1)^(i_k + j_k) \
        &= cases(
            1 " if " i_k = 0\, j_k = 0,
            0 " if " i_k = 0 \, j_k =1,
            0 " if " i_k = 1 \, j_k =0,
            1 " if " i_k = 1\, j_k = 1,
        ) \
        &= delta_(i_k , j_k).
    $
    Then we have that the total product is
    $
        EE_(Lambda_G) bra(i) Lambda_G ket(i) bra(j) Lambda_G ket(j) = product_(k = 1)^n delta_(i_k , j_k) = delta_(i,j).
    $
]

= Haar Integrals <sec_tsp_appendix>

In this section we present the more technical work needed to state our results in @sec_tsp_weak_coupling. @lem_two_heisenberg_interactions and @lem_sandwiched_interaction are used to compute the effects of the randomized interactions in a form that are usable in the main result of @lem_tsp_transitions. @lem_haar_two_moment can be derived from Appendix C in @brandao2021complexity.

#lemma()[
    Let $EE_U$ denote the expectation over the Haar measure over the set of unitary matrices acting on a $dim$ dimensional Hilbert space. Then for $ket(i_1), ket(i_2), ..., ket(k_2)$ drawn from an orthonormal basis
    $
        &EE_U [ bra(i_1) U ket(j_1) bra(i_2) U ket(j_2) bra(k_1) U^dagger ket(l_1) bra(k_12) U^dagger ket(l_2) ] \
        =& 1 / (dim^2 - 1) (delta_(i_1, l_1) delta_(j_1, k_1) delta_(i_2, l_2) delta_(j_2, k_2) + delta_(i_1, l_2) delta_(j_1, k_2) delta_(i_2, l_1) delta_(j_2, k_1)) \
        &- 1 / (dim(dim^2 - 1)) (delta_(i_1, l_2) delta_(j_1, k_1) delta_(i_2, l_1) delta_(j_2, k_2) + delta_(i_1, l_1) delta_(j_1, k_2) delta_(i_2, l_2) delta_(j_2, k_1)).
    $ <eq_haar_two_moment_integral>
] <lem_haar_two_moment>

#lemma()[
    Let $G(t)$ denote the Heisenberg evolved random interaction $G(t) = e^(i H t) G e^(-i H t)$ for a total Hamiltonian $H$. After averaging over the interaction measure the product $G(x) G(y)$ can be computed as
    $
        EE_G [G(x) G(y)] = 1 / (dim + 1) (sum_((i,j), (k,l)) e^(Delta (i, j | k, l) (x - y)) ketbra(i\, j, i\, j) + id ).
    $
] <lem_two_heisenberg_interactions>
#proof()[
    The overall structure of this proof is to evaluate the product in the Hamiltonian eigenbasis and split the product into three factors: a phase contribution from the time evolution, a Haar expectation from the eigenvectors of the random interaction, and the eigenvalue expectation of the random interaction. Since this involves the use of multiple indices, it will greatly simplify the proof to use a single index over the total Hilbert space $hilb$ as opposed to two indices over $hilb_S tp hilb_E$. For example, the index $a$ should be thought of as a pair $(a_s, a_e)$, and functions $lambda (a)$ should be thought of as $lambda (a_s, a_e)$. Once the final form of the expression is reached we will substitute in pairs of indices for easier use of the lemma in other places.
    $
        EE_G [ G(x) G(y) ] &= EE_(Lambda_G) EE_(U_G) e^(+i H x) U_G Lambda_G U_G^dagger e^(-i H x) e^(+i H y) U_G Lambda_G U_G^dagger e^(-i H y) \
        &= EE_(Lambda_G) EE_(U_G) [sum_a e^(+i lambda(a)x) ket(a) bra(a) U_G sum_b Lambda_G (b) ket(b) bra(b) U_G^dagger \
            &quad sum_c e^(-i lambda(c) (x - y)) ket(c) bra(c) U_G sum_d Lambda_G (d) ket(d) bra(d) U_G^dagger sum_e e^(-i lambda(e) y) ket(e) bra(e) ] \
        &= sum_(a,b,c,d,e) ket(a) bra(e) e^(-i (lambda(c) - lambda(a))x) e^(-i (lambda(e) - lambda(c))y) \
        &quad times EE_(U_G) [ bra(a) U_G ket(b) bra(c) U_G ket(d) bra(b) U_G^(dagger) ket(c) bra(d) U_G^dagger ket(e) ] EE_(Lambda_G) [ Lambda_G (b) Lambda_G (d) ] \
        &= sum_(a, b, c, d, e) delta_(b d) ket(a) bra(e) e^(-i (lambda(c) - lambda(a))x) e^(-i (lambda(e) - lambda(c))y) \
        &quad times EE_(U_G) [ bra(a) U_G ket(b) bra(c) U_G ket(d) bra(b) U_G^(dagger) ket(c) bra(d) U_G^dagger ket(e) ].
    $ <eq_appendix_tmp_1>
    We used the fact that the eigenvalues of $G$ are I.I.D with variance 1 to make the substitution $EE_(Lambda_G) [Lambda_G (b) Lambda_G (d)] = delta_(b, d)$. This allows us to reduce the sum over $d$ to the condition when $d=b$, which greatly simplifies the Haar expectation we have to take. As our eigenvectors are chosen from the Clifford group, which forms a 2-design, we can use @lem_haar_two_moment to write
    $
        EE_(U_G) [bra(a) U_G ket(b) bra(c) U_G ket(b) bra(b) U_G^(dagger) ket(c) bra(b) U_G^dagger ket(e) ] &= EE_U [bra(a) U ket(b) bra(c) U ket(b) bra(b) U^(dagger) ket(c) bra(b) U^dagger ket(e) ] \
        &=1 / (dim^2 - 1) (delta_(a c) delta_(c e) + delta_(a e) - frac(1, "dim") (delta_(a c) delta_(c e) + delta_(a e))). "      "
    $
    Plugging this into @eq_appendix_tmp_1 above yields the following
    $
        &= frac(1, "dim"^2 - 1) sum_(a, b, c, e) ket(a) bra(e) e^(-i (lambda(c) - lambda(a))x) e^(-i (lambda(e) - lambda(c))y) (delta_(a c) delta_(c e) + delta_(a e) - frac(1, "dim") (delta_(a c) delta_(c e) + delta_(a e))) \
        &= frac(1, "dim"^2 - 1) (1 - frac(1, "dim")) sum_(a, b, c, e) ket(a) bra(e) e^(-i (lambda(c) - lambda(a))x) e^(-i (lambda(e) - lambda(c))y) delta_(a e) (1 + delta_(a c)) \
        &= frac(1, "dim"^2 - 1) (1 - frac(1, "dim")) sum_(a, b, c) ket(a) bra(a) e^(i (lambda(a) - lambda(c))(x-y)) (1 + delta_(a c)) \
        &= frac(1 ("dim" - 1), "dim"^2 - 1) sum_(a,c) ket(a) bra(a) e^(i (lambda(a) - lambda(c))(x - y)) (1 + delta_(a c)) \
        &= frac(1, "dim" + 1) (sum_(a,c) e^(i (lambda(a) - lambda(c))(x-y)) ket(a) bra(a) + identity).
    $

    Reindexing by $a |-> i,j$, $c |-> k,l$, and plugging in the definition of $Delta$ yields the statement of the lemma.
]

#lemma()[
    Given two Heisenberg evolved random interactions $G(x)$ and $G(y)$ we can compute their action on the outer product $ketbra(i\, j, k\, l)$ as
    $
        &EE_(G) [ G(x) ketbra(i\, j, k\, l) G(y)] \
        &= 1 / (dim + 1) (ketbra(i\, j, k\, l) + braket(i\, j, k\, l) sum_(m,n) e^(Delta (m,n | i,j) (x - y)) ketbra(m\, n, m\, n)).
    $
] <lem_sandwiched_interaction>
#proof()[
    This proof is structured the same as @lem_two_heisenberg_interactions and similarly we will use a single index of the total Hilbert space $hilb$ and switch to two indices to match the rest of the exposition.

    $
        EE_G [ G(x) ket(a) bra(b) G(y) ] &= EE_G [ e^(i H x) U_G D U_G^(dagger) e^(-i H x) ket(a) bra(b) e^(i H y) U_G D U_G^dagger e^(-i H y) ] \
        &= sum_(c, d, e, f) e^(i (lambda(c) - lambda(a))x) e^(i (lambda(b) - lambda(f))y) \
        &quad times EE_G [ ket(c) bra(c) U_G D(d) ket(d) bra(d) U_G^dagger ket(a) bra(b) U_G D(e) ket(e) bra(e) U_G^dagger ket(f) bra(f) ] \
        &= sum_(c, d, e, f) e^(i (lambda(c) - lambda(a))x) e^(i (lambda(b) - lambda(f))y) ket(c) bra(f) \
        &quad times EE_(Lambda_G) [Lambda_G (d) Lambda_G (e) ] EE_(U_G) [ bra(c) U_G ket(d) bra(b) U_G ket(e) bra(d) U_G^dagger ket(a) bra(e) U_G^dagger ket(f) ] \
        &= sum_(c,d,f) e^(i (lambda(c) - lambda(a))x) e^(i (lambda(b) - lambda(f))y) ket(c) bra(f) \
        &quad times EE_(U_G) [ bra(c) U_G ket(d) bra(b) U_G ket(d) bra(a) overline(U_G) ket(d) bra(f) overline(U_G) ket(d) ] \
        &= frac(1, "dim"^2 - 1) sum_(c,d,f) e^(i (lambda(c) - lambda(a))x) e^(i (lambda(b) - lambda(f))y) ket(c) bra(f) (delta_(c a) delta_(b f) + delta_(c f) delta_(a b))(1 - frac(1, "dim")) \
        &= frac(1, "dim" + 1) sum_(c,f) e^(i (lambda(c) - lambda(a))x) e^(i (lambda(b) - lambda(f))y) ket(c) bra(f) (delta_(c a) delta_(b f) + delta_(c f) delta_(a b)) \
        &= frac(1, "dim" + 1) (ket(a) bra(b) + delta_(a b) sum_(c) e^(i(lambda(c) - lambda(a))(x-y)) ket(c) bra(c)).
    $
    We used the fact that $EE_(Lambda_G) [Lambda_G (d) Lambda_G (e)] = delta_(d, e)$ to eliminate the sum over $e$. We used the same Haar 2-design argument to compute the expectation value $EE_(U_G) [dot]$. Re-indexing by $a |-> (i,j)$, $b |-> (k,l)$ and $c |-> (m,n)$ results in the expression given in the statement of the lemma.
]

Now that we have computed the expected Heisenberg evolution of a two-body interaction term we can compute the second order transition amplitudes. We will not restate the lemma here for brevity.
#proof([of @lem_tsp_transitions])[
    To start we would like to note that we will use a single index notation to refer to the joint system-environment eigenbasis during this proof to help shorten the already lengthy expressions. We will convert back to a double index notation to match the statement of the theorem. We start from the expression for the first derivative of the channel $frac(diff, partial alpha) Phi_G (rho_S)$ given by @eq_first_order_alpha_derivative. To take the second derivative there are six factors involving $alpha$, so we will end up with six terms. We repeat @eq_first_order_alpha_derivative below, add a derivative, and label each factor containing an $alpha$ for easier computation
    $
        frac(partial^2, partial alpha^2) Phi_G (rho_S) &= frac(partial, partial alpha) (integral_(0)^(1) underbrace(e^(i s (H+alpha G)t), "(A)") (i t G) underbrace(e^(i (1-s) (H+alpha G)t), "(B)") d s " " rho " " underbrace(e^(-i(H+alpha G)t), "(C)") ) \
        &quad +frac(partial, partial alpha) ( underbrace(e^(i(H+alpha G)t), "(D)") " " rho " " integral_(0)^1 underbrace(e^(-i s (H+alpha G) t), "(E)") (- i t G) underbrace(e^(-i (1-s) (H+alpha G)t), "(F)") d s ).
    $ <eq_second_derivative_labels>

    Our goal is to get each of these terms in a form in which we can use either @lem_two_heisenberg_interactions or @lem_sandwiched_interaction.
    $
        (A) &= i t integral_0^1 (frac(partial, partial alpha) e^(i s_1 (H+ alpha G)t)) G e^(i(1-s_1)(H+alpha G)t) d s_1 rho e^(-i (H+alpha G)t) |_(alpha=0) \
        &= (i t)^2 integral_0^1 (integral_0^1 e^(i s_1 s_2 (H+alpha G)t) s_1 G e^(i s_1 (1-s_2) (H+alpha G)t) d s_2) G e^(i(1-s_1) (H+alpha G)t) d s_1 rho e^(-i(H+alpha G) t) |_(alpha=0) \
        &= -t^2 integral_0^1 integral_0^1 e^(i s_1 s_2 H t) G e^(-i s_1 s_2 H t) e^(i s_1 H t) G e^(-i s_1 H t) s_1 d s_1 d s_2 e^(i H t) rho e^(-i H t) \
        &= -t^2 integral_0^1 integral_0^1 G(s_1 s_2 t) G(s_1 t) s_1 d s_1 d s_2 rho(t).
    $ <eq_second_deriv_alpha_first_term>

    $
        (B) &= i t integral_0^1 e^(i s_1 (H + alpha G)t) G frac(partial, partial alpha)(e^(i(1-s_1)(H + alpha G)t)) d s_1 rho e^(-i(H + alpha G) t) |_(alpha = 0) \
        &= (i t)^2 integral_0^1 e^(i s_1 (H + alpha G)t) G (integral_0^1 e^(i(1-s_1)s_2 (H + alpha G)t) (1-s_1) G e^(i(1 - s_1)(1 - s_2)(H + alpha G)t) d s_2) d s_1 " " rho e^(-i ( H + alpha G)t) |_(alpha = 0) \
        &= -t^2 integral_0^1 integral_0^1 e^(i s_1 H t) G e^(i(1-s_1)s_2 H t) G e^(i(1-s_1)(1-s_2) H t) (1-s_1) d s_1 d s_2 " " rho e^(-i H t) \
        &= -t^2 integral_0^1 integral_0^1 e^(i s_1 H t) G e^(-i s_1 H t) e^(i(s_1 + s_2 - s_1 s_2) H t) G e^(-i (s_1 + s_2 - s_1 s_2) H t) (1-s_1) d s_1 d s_2 " " rho(t) \
        &= -t^2 integral_0^1 integral_0^1 G(s_1 t) G((s_1 + s_2 - s_1 s_2)t) (1-s_1) d s_1 d s_2 " " rho(t)
    $

    $
        (C) &= i t integral_0^1 e^(i s (H + alpha G)t) G e^(i(1-s) (H + alpha G) t) d s " "rho " " frac(partial, partial alpha) ( e^(-i (H + alpha G) t) ) |_(alpha = 0) \
        &= (i t) (-i t) integral_0^1 e^(i s (H + alpha G)t) G e^(i (1 - s) (H + alpha G)t) d s " " rho " " ( integral_0^1 e^(-i s (H + alpha G)t) G e^(-i (1- s) ( H + alpha G)t ) d s)|_(alpha = 0) \
        &= + t^2 (integral_0^1 e^(i s H t) G e^(-i s H t) d s) e^(i H t) rho e^(-i H t) (integral_0^1 e^(i (1-s) H t) G e^(-i (1-s) H t) d s) \
        &= + t^2 integral_0^1 G(s t) d s " " rho(t) integral_0^1 G((1-s)t) d s
    $

    $
        (D) &= (-i t) frac(partial, partial alpha) (e^(i(H + alpha G)t)) rho integral_0^1 e^(-i s (H + alpha G)t) G e^(-i (1-s)(H + alpha G)t) d s |_(alpha = 0) \
        &= t^2 (integral_0^1 e^(i s (H+ alpha G)t) G e^(i (1-s) (H + alpha G)t)d s) rho integral_0^1 e^(-i s (H + alpha G)t) G e^(-i (1-s)(H + alpha G)t) d s |_(alpha = 0) \
        &= t^2 integral_0^1 e^(i s H t) G e^(-i s H t) d s " " rho (t) integral_0^1 e^(i (1-s) H t) G e^(-i (1-s) H t) d s \
        &= t^2 integral_0^1 G(s t) d s " " rho(t) " " integral_0^1 G((1-s)t) d s
    $

    $
        (E) &= (-i t) e^(i (H+ alpha G) t) " " rho " " integral_0^1 frac(partial, partial alpha) (e^(-i s_1 (H + alpha G)t)) G e^(-i (1-s_1)(H + alpha G)t) d s_1 |_(alpha = 0) \
        &= - t^2 e^(i(H + alpha G)t) " " rho " " integral_0^1 (integral_0^1 e^(-i s_1 s_2 (H + alpha G) t) (s_1 G) e^(-i s_1 (1-s_2) (H + alpha G)t) d s_2) G e^(-i(1-s_1)(H + alpha G)t) d s_1 |_(alpha = 0) \
        &= -t^2 e^(i H t) rho e^(-i H t) integral_0^1 integral_0^1 e^(i (1 - s_1 s_2) H t) G e^(-i (s_1 - s_1 s_2)H t) G e^(-i (1-s_1)H t) s_1 d s_1 d s_2 \
        &= -t^2 rho(t) integral_0^1 integral_0^1 G((1- s_1 s_2) t) G((1-s_1)t) s_1 d s_1 d s_2
    $

    $
        (F) &= (-i t) e^(i(H + alpha G) t) rho integral_0^1 e^(-i s_1 ( H + alpha G) t) G frac(partial, partial alpha) ( e^(-i (1-s_1) ( H +alpha G)t)) d s_1 |_(alpha = 0) \
        &= (-i t)^2 e^(i (H + alpha G)t) rho integral_0^1 e^(-i s_1 (H + alpha G)t) G (integral_0^1 e^(-i(1-s_1) s_2 (H + alpha G)t) (1-s_1) G e^(-i(1-s_1) (1-s_2) (H + alpha G) t) d s_2) d s_1 |_(alpha = 0) \
        &= -t^2 e^(-i H t) rho e^(-i H t) integral_0^1 integral_0^1 e^(i (1- s_1) H t) G e^(-i (1-s_1) H t) e^(i (1-s_1)(1-s_2) H t) G e^(-i(1-s_1)(1-s_2) H t) (1-s_1) d s_1 d s_2 \
        &= -t^2 rho(t) integral_0^1 integral_0^1 G((1-s_1)t) G((1-s_1)(1 - s_2) t) (1-s_1)d s_1 d s_2
    $

    Now our goal is to compute the effects of averaging over the interaction $G$ on the above terms, starting with $(A)$. As this involves a lot of index manipulations, similarly to the proofs of Lemmas @lem_two_heisenberg_interactions and @lem_sandwiched_interaction we will use a single index for the total system-environment Hilbert space and switch back to a double index to state the results. We will make heavy use of Lemma @lem_two_heisenberg_interactions.
    $
        EE_G (A) &= -t^2 integral_0^1 integral_0^1 EE_G [ G(s_1 s_2 t) G(s_1 t) ] s_1 d s_1 d s_2 rho(t) \
        &= frac(-t^2, dim + 1) integral_0^1 integral_0^1 (sum_(i,j) e^(i (lambda(i) - lambda(j)) (s_1 s_2 t - s_1 t)) ket(i)bra(i) + identity) s_1 d s_1 d s_2 rho(t) \
        &= frac(- t^2, dim + 1) (sum_(i) sum_(j : lambda(i) != lambda(j)) integral_0^1 integral_0^1 e^(i(lambda(i) - lambda(j))t (s_1 s_2 - s_1)) s_1 d s_1 d s_2 ket(i)bra(i) + sum_(i) sum_(j : lambda(i) = lambda(j))frac(1, 2) ket(i)bra(i) + frac(1, 2) identity) rho(t) \
        &= frac(- t^2, dim + 1) (sum_i sum_(j : lambda(i) != lambda(j)) frac(1 - i (lambda(i) - lambda(j))t - e^(-i (lambda(i) - lambda(j))t), t^2 (lambda(i) - lambda(j))^2) ket(i)bra(i) + frac(1, 2) sum_(i) (eta(i) + 1) ket(i)bra(i) ) rho(t) \
        &= frac(- 1, dim + 1)(sum_(i) sum_(j: Delta_(i j) != 0) frac(1 - i Delta_(i j)t - e^(-i Delta_(i j) t), Delta_(i j)^2) ket(i)bra(i) + frac(t^2, 2) sum_(i) (eta(i) + 1)ket(i)bra(i) ) rho(t)
    $

    We can similarly compute the averaged $(B)$ term:
    $
        EE_G (B) &= -t^2 integral_0^1 integral_0^1 EE_G [ G(s_1 t) G((s_1 + s_2 - s_1 s_2) t) ] (1-s_1) d s_1 d s_2 " " rho(t) \
        &= frac(- t^2, dim + 1) integral_0^1 integral_0^1 (sum_(i,j) e^(i (lambda(i) - lambda(j))(s_1 s_2 - s_2) t) ket(i)bra(i) + identity) (1 -s_1) d s_1 d s_2 " " rho \
        &= frac(- t^2, dim + 1) (sum_(i) sum_(j : lambda(i) != lambda(j)) integral_0^1 integral_0^1 e^(i(lambda(i) - lambda(j))t (s_1 s_2 - s_2)} (1 - s_1) d s_1 d s_2 ket(i)bra(i) + sum_(i) sum_(j : lambda(i) = lambda(j))frac(1, 2) ket(i)bra(i) + frac(1, 2) identity) rho(t) \
        &= frac(- t^2, dim + 1) (sum_i sum_(j : lambda(i) != lambda(j)) frac(1 - i (lambda(i) - lambda(j))t - e^(-i (lambda(i) - lambda(j))t), t^2 (lambda(i) - lambda(j))^2) ket(i)bra(i) + frac(1, 2) sum_(i) (eta(i) + 1) ket(i)bra(i) ) rho(t) \
        &= frac(-1, dim + 1)(sum_(i) sum_(j: Delta_(i j) != 0) frac(1 - i Delta_(i j)t - e^(-i Delta_(i j) t), Delta_(i j)^2) ket(i)bra(i) + frac(t^2, 2) sum_(i) (eta(i) + 1)ket(i)bra(i) ) rho(t),
    $
    which we note is identical to $EE_G (A)$. As terms $(C)$ and $(D)$ involve a different method of computation we skip them for now and compute $(E)$ and $(F)$.
    $
        EE_G (E) &= -t^2 rho(t) integral_0^1 integral_0^1 EE_G [ G((1- s_1 s_2) t) G((1-s_1)t) ] s_1 d s_1 d s_2 \
        &= frac(- t^2, dim + 1) rho(t) integral_0^1 integral_0^1 (sum_(i,j) e^(i(lambda(i) - lambda(j)) t (s_1 - s_1 s_2)) ket(i)bra(i) + identity ) s_1 d s_1 d s_2 \
        &= frac(- t^2, dim + 1) rho(t) (sum_i sum_(j : lambda(i) != lambda(j)) frac(1 + i (lambda(i) - lambda(j))t - e^(i(lambda(i) - lambda(j))t), t^2 (lambda(i) - lambda(j))^2)ket(i)bra(i) + frac(1, 2) sum_(i) (eta(i) + 1 )ket(i)bra(i)) \
        &= frac(- 1, dim + 1) rho(t) (sum_i sum_(j: (Delta_(i j) != 0)) frac(1 + i Delta_(i j)t - e^(i Delta_(i j)t), Delta_(i j)^2) ket(i)bra(i) + frac(t^2, 2)sum_i (eta(i) + 1) ket(i)bra(i)).
    $
    Computing $(F)$ yields
    $
        EE_G (F) &= -t^2 rho(t) integral_0^1 integral_0^1 EE_G [ G((1-s_1)t) G((1-s_1)(1 - s_2) t) ] (1-s_1)d s_1 d s_2 \
        &= frac(- t^2 sigma^2, dim + 1) rho(t) integral_0^1 integral_0^1 (sum_(i,j) e^(i(lambda(i) - lambda(j))t (s_2 - s_1 s_2))ket(i)bra(i) + identity) (1-s_1) d s_1 d s_2 \
        &= frac(- t^2, dim + 1) rho(t) (sum_(i) sum_(j : lambda(i) != lambda(j)) frac(1 + i (lambda(i) - lambda(j))t - e^(i (lambda(i) - lambda(j))t), t^2 (lambda(i) - lambda(j))^2) ket(i)bra(i) +frac(1, 2) sum_(i) (eta(i) + 1) ket(i)bra(i)) \
        &= frac(- 1, dim + 1) rho(t) (sum_i sum_(j: (Delta_(i j) != 0)) frac(1 + i Delta_(i j)t - e^(i Delta_(i j)t), Delta_(i j)^2) ket(i)bra(i) + frac(t^2, 2)sum_i (eta(i) + 1) ket(i)bra(i))
    $
    which is identical to $EE_G (E)$.

    The last two terms $(C) = (D)$ are computed as follows:
    #set math.equation(number-align: top)
    $
        EE_G (C) &= t^2 integral_0^1 integral_0^1 EE_G [ G(s_1 t) rho(t) G((1-s_2)t) ] " " d s_1 d s_2 \
        &= t^2 sum_(i,j) rho_(i j) e^(i(lambda(i) - lambda(j))t) integral_0^1 integral_0^1 EE_G [ G(s_1 t) ket(i)bra(j) G((1-s_2)t) ] " " d s_1 d s_2 \
        &= frac(t^2, dim + 1) sum_(i,j) rho_(i j) e^(i(lambda(i) - lambda(j))t) ( ket(i)bra(j) + delta_(i j) sum_(a) integral_0^1 integral_0^1 e^(i(lambda(a) - lambda(i))(s_1 + s_2 - 1)t) d s_1 d s_2 ket(a)bra(a)) \
        &= frac(t^2, dim + 1) sum_(i,j) rho_(i j) e^(i Delta_(i j) t) (ket(i)bra(j) + delta_(i j) sum_(a : Delta_(a i) != 0) frac(2(1- cos (Delta_(a i) t)), Delta_(a i)^2 t^2) ket(a)bra(a) + delta_(i j) sum_(a : Delta_(a i) = 0) ket(a)bra(a))
    $


    We can now combine each of these terms to offer the full picture of the output of the channel to second order. We make two modifications to the results from each sum: first, we will switch to double index notation to make for easier use in other areas, and secondly we let $rho = ket(i\, j)bra(k\, l)$. We note that the first term in the following equation is provided by $(A) + (B)$, the second through $(E) + (F)$, and the last two through $(C) + (D)$.
    $
        &EE_G [ frac(partial^2, partial alpha^2) Phi_G(ket(i\, j)bra(k\, l))|_(alpha = 0) ] \
        &= -frac(2 e^(i Delta(i, j|k, l) t), dim + 1) (sum_((a,b): Delta(i, j|a, b) != 0) frac(1 - i Delta(i, j|a, b)t - e^(-i Delta(i, j|a, b) t), Delta(i, j|a, b)^2) \
            &+ sum_((a,b): Delta(k, l|a, b) != 0) frac(1 + i Delta(k, l|a, b) t - e^(i Delta(k, l|a, b) t), Delta(k, l|a, b)^2) + frac(t^2, 2)(eta(i, j) + eta(k, l)) ) ket(i\,j)bra(k\,l) \
        & +delta_(i,k) delta_(j,l) frac(2 e^(i Delta(i, j|k, l)t), dim+1) ( sum_((a,b): Delta(i, j|a, b) != 0 ) frac(2(1- cos (Delta(i, j|a, b)t)), Delta(i, j|a, b)^2) ket(a\,b)bra(a\,b) + t^2 sum_((a,b) : Delta(i, j|a, b) = 0) ket(a\,b)bra(a\,b))
    $ <eq_second_order_output>
    #set math.equation(number-align: bottom)

    The last step we need is to use the half angle formula to change the cosine to a sine
    $
        frac(2(1 - cos(Delta(i, j| a, b)t)), Delta(i, j|a, b)^2) = frac(2(1 - (1 - 2 sin^2(frac(Delta(i, j|a, b)t, 2)))), Delta(i, j|a, b)^2) = t^2 sinc^2 (frac(Delta(i, j|a, b) t, 2)), "     "
    $ <eq_trig_end>
    which yields the statement.

    We can compute these by plugging in to Eq. @eq_el_gigante again, which yields
    $
        &EE_G [ bra(i'\, j') cal(T) ( ket(i\, j)bra(i\, j) ) ket(i'\, j') ] = cases(
            tilde(alpha)^2 sinc^2(Delta(i, j | i', j') t / 2) & (i, j) != (i', j') \
            -tilde(alpha)^2 sum_((a,b) != (i, j)) sinc^2(Delta(a, b|i, j) t / 2) & (i,j) = (i', j')
        ). "  "
    $ <eq_system_environment_transitions>

    The $(i, j) != (i', j')$ case should be apparent, the first term with the coherence factors $chi$ are zero and the second term is what remains. The $(i,j) = (i', j')$ case can be seen as follows. For the first term we have
    $
        - frac(alpha^2 e^(i Delta(i, j| i, j) t), dim + 1)(chi(i, j) + chi(i, j)^* + frac(t^2, 2)(eta(i, j) + eta(i, j))) ket(i\,j)bra(i\,j).
    $
    We first compute the sum $chi (i,j) + chi (i,j)^*$ as
    $
        chi(i, j) + chi(i, j)^* &= sum_(a,b: Delta(i, j, |a, b) != 0) frac(1 - i Delta(i, j|a, b)t - e^(-i Delta(i, j|a, b) t), Delta(i, j|a, b)^2) \
        &quad+ sum_(a,b: Delta(i, j, |a, b) != 0) frac(1 + i Delta(i, j|a, b)t - e^(+i Delta(i, j|a, b) t), Delta(i, j|a, b)^2) \
        &= sum_(a,b: Delta(i, j| a, b) != 0) frac(2 - e^(-i Delta(i, j| a, b) t) - e^(+i Delta(i, j| a, b) t), Delta(i, j|a, b)^2) \
        &= sum_(a,b: Delta(i, j| a, b) != 0) t^2 sinc^2 ( frac(Delta(i, j| a, b) t, 2) ),
    $
    where the last step follows from a trigonometric identity (see @eq_trig_end). Since $sinc(0) = 1$ the $eta(i, j)$ term can be expressed as $eta(i, j) = sum_(a,b : Delta(i, j|a, b) = 0) sinc^2 ( frac(Delta(i, j| a, b) t, 2) )$. Plugging this into Eq. @eq_el_gigante gives
    $
        &EE_G [ bra(i\,j) cal(T) (ket(i\,j)bra(i\,j)) ket(i\,j) ] \
        &= bra(i\,j) (-frac(alpha^2 t^2, dim + 1) sum_(a,b) sinc^2 ( frac(Delta(i, j| a, b) t, 2) ) ket(i\,j)bra(i\,j) + sum_(a,b) sinc^2( frac(Delta(i, j | a, b)t, 2) ) ket(a\,b)bra(a\,b) ) ket(i\,j) \
        &= -frac(alpha^2 t^2, dim + 1) sum_((a,b) != (i,j)) sinc^2 ( frac(Delta(i, j| a, b) t, 2) ).
    $
    As a by-product of this computation we have also shown that $tr(cal(T)(rho)) = 0$ and that our mapping is trace preserving to $O(alpha^2)$.
]

#proof([of @thm_remainder_bound])[
    First we note that although $R_Phi (rho) = alpha^3 / 6 diff_alpha^3 Phi(rho) |_(alpha = alpha_star)$ for a specific value $alpha_star > 0$ our proof will actually hold for any value of $alpha_star > 0$. To compute the trace norm we will use the triangle inequality, unitary invariance of the Schätten norms, and submultiplicativity. To start,
    $
        norm(diff_alpha^3 Phi(rho))_1 &= norm((diff^3) / (diff alpha^3) EE_G tr_E e^(i (H + alpha G)t) rho tp rho_E e^(-i (H + alpha G)t))_1 \
        &<= EE_G norm((diff^3) / (diff alpha^3) e^(i (H + alpha G)t) rho tp rho_E e^(-i (H + alpha G)t))_1,
    $
    where we can take $EE_G$ out of the norm via the triangle inequality and we can remove the trace via Proposition 1 of @rastegin2012relations, which proves $norm(tr_E [X])_1 <= norm(X)_(dim_E) <= norm(X)_1$. To proceed we use the decomposition of the second derivatives from the proof of @lem_tsp_transitions, specifically @eq_second_derivative_labels. This gives the following
    $
        norm(R_Phi)_1 &<= alpha^3 / 6 EE_G norm(diff_alpha ((A) + (B) + (C) + (D) + (E) + (F))|_(alpha = alpha_star))_1 \
        &<= alpha^3 / 6 (EE_G norm(diff_alpha (A)|_(alpha = alpha_star))_1 + ... + EE_G norm(diff_alpha (F)|_(alpha = alpha_star))_1).
    $
    We will demonstrate how this can be computed for the first term $diff_alpha (A)$. Using @eq_second_deriv_alpha_first_term and letting $H_alpha = H + alpha G$ for brevity we can write
    $
        diff_alpha (A) = -t^2 diff_alpha integral_0^1 integral_0^1 underbrace(e^(i s_1 s_2 H_alpha t), "(A.1)") G underbrace(e^(i s_1(1-s_2) H_alpha t), "(A.2)") G underbrace(e^(i(1 - s_1) H_alpha t), "(A.3)") rho underbrace(e^(-i H_alpha t), "(A.4)") s_1 d s_1 d s_2,
    $
    where there are four spots for the derivative to act via Duhamel's formula. We will show only one of these terms, starting with (A.1)
    #set math.equation(number-align: horizon)
    $
        &"(A.1)" = -t^2 integral_0^1 integral_0^1 diff_alpha (e^(i s_1 s_2 H_alpha t)) G e^(i s_1(1-s_2) H_alpha t) G e^(i(1 - s_1) H_alpha t) rho e^(-i H_alpha t) s_1 d s_1 d s_2 \
        &= (i t)^3 integral_0^1 integral_0^1 integral_0^1 e^(i s_1 s_2 s_3 H_alpha t) G e^(i s_1 s_2 (1 - s_3) H_alpha t) G e^(i s_1 (1 - s_2) H_alpha t) G e^(i (1-s_1) H_alpha t ) rho e^(-i H_alpha t) s_1^2 s_2 d s_1 d s_2 d s_3.
    $
    #set math.equation(number-align: bottom)
    Our goal is to compute the 1-norm of the above expression at $alpha = alpha_star$. We can do so using the triangle inequality to move the norms into the integrand and then use submultiplicativity and unitary invariance to achieve
    $
        &norm(e^(i s_1 s_2 s_3 H_alpha_star t) G e^(i s_1 s_2 (1 - s_3) H_alpha_star t) G e^(i s_1 (1 - s_2) H_alpha_star t) G e^(i (1-s_1) H_alpha_star t ) rho e^(-i H_alpha_star t))_1 \
        &<= norm(e^(i s_1 s_2 s_3 H_alpha_star t) G)_1 norm(e^(i s_1 s_2 (1 - s_3) H_alpha_star t) G)_1 norm(e^(i s_1 (1 - s_2) H_alpha_star t) G)_1 norm(e^(i (1-s_1) H_alpha_star t ) rho e^(-i H_alpha_star t))_1 \
        &<= norm(G)_1^3 norm(rho)_1 = norm(G)_1^3.
    $
    Similar computations can be carried out for the other three terms (A.2) - (A.4). In total these yield the inequality
    $
        alpha^3 / 6 integral norm(diff_alpha (A))_1 d G &<= (alpha t)^3 / 6 EE_G integral_0^1 integral_0^1 integral_0^1 norm(G)_1^3 (s_1^2 s_2 + s_1^2(1 - s_2) + s_1(2 - s_1)) d s_1 d s_2 d s_3 \
        &<= 4 / 6 (alpha t)^3 EE_G norm(G)_1^3.
    $ <eq_remainder_bound_on_A>

    #h(5mm) Now that we have computed a bound for the norm of the derivative acting on $(A)$ we only have terms $(B)$ through $(F)$ to compute. These can all be checked to satisfy the same bound on $(A)$ from @eq_remainder_bound_on_A, and as there are six terms in total we have the inequality
    $
        norm(R_Phi)_1 <= 4 (alpha t)^3 EE_G norm(G)_1^3,
    $
    which holds for all inputs $rho$.

    Our last remaining problem is to compute the expected norm of $G$. Using the decomposition of our interaction $G = U_G Lambda_G U_G^dagger$ to get
    $
        EE_G norm(G)_1^3 = EE_(Lambda_G) EE_(U_G) norm(U_G Lambda_G U_G^dagger)_1^3 = EE_(Lambda_G) norm(Lambda_G)_1^3 = EE_(Lambda_G) sum_(i = 1)^(dim) abs(Lambda_G (i))^3 = dim ,
    $
    Since $Lambda_G$ is just $plus.minus 1$ times a Pauli $Z$ string each eigenvalue has norm 1. This gives the final inequality
    $
        norm(R_Phi)_1 <= 4 dim (alpha t)^3.
    $
]

// = Scratch

// == Template Thermal State Prep Proof
// I'm thinking of including a "template" theorem that can be used to simplify the 4 proofs contained in the following section. Let $rho_"fix"$ denote the unique fixed point for a channel $EE_gamma [id + cal(T)_"on"^((gamma))]$ and $tilde(lambda_star)$ the spectral gap of the scaled transition matrix, so $tilde(alpha)^2 tilde(lambda_star) = lambda_star$. Then we have
// $
//     norm(rho_"fix" - (EE_gamma Phi_gamma)^(compose L) (rho))_1 &<= norm(rho_"fix" - (id + EE_gamma cal(T)_"on"^((gamma)))^(compose L))_1 + L norm(EE_gamma cal(T)_"off"^((gamma)) + R_Phi)_1 \
//     &<= norm(rho_"fix" - (id + EE_gamma cal(T)_"on"^((gamma)))^(compose L))_1 + L (norm(EE_gamma cal(T)_"off"^((gamma)))_1 + norm(R_Phi)_1) \
//     &<= norm(rho_"fix" - (id + EE_gamma cal(T)_"on"^((gamma)))^(compose L))_1 + L ((8 alpha^2) / delta_min^2 + (16 sqrt(2)) / sqrt(pi) dim_S (alpha t)^3).
// $
// So now in order to balance these terms we can set $alpha = 1\/(dim_S delta_min^2 t^3 )$ and the expression on the right becomes $L alpha^2 / delta_min^2 (8 + 16 sqrt(2 / pi)).$ Now using Jerison's theorem we can argue that
// $
//     L >= dim^2 / (alpha^2 t^2 tilde(lambda_star)) J
// $
// is sufficient to guarantee that the distance to the fixed point is $tilde(O)(epsilon)$.
// Now we note that the right hand side forces us to require $L alpha^2 / delta_min^2 in tilde(O)(epsilon)$ holds only if
// $
//     (dim^2) / (delta_min^2 t^2 tilde(lambda_star) ) in tilde(O)(epsilon)
// $
// can be satisfied if $t = dim / (delta_min sqrt(epsilon tilde(lambda_star)))$. and then we are done.

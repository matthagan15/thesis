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

#heading("Preparing Thermal Quantum States", level: 1, supplement: "Chapter") <ch:thermal_state_prep>
Thermal states of the form $e^(-beta H) / cal(Z)$ are ubiquitous in physics. These are the states that we believe physical systems take whenever they are cooled (or heated) to an inverse temperature of $beta$. They are typically used to estimate observables $O$ of interest in physically relevant states, oftentimes of the form $angle.l O angle.r_beta = tr(O e^(-beta H)/ cal(Z))$. These observables could be anything from dipole moments, magnetizations, or two body correlators. Classically these states correspond to the canonical ensemble, a distribution over phase space that tells us the probability of finding a particle at a particular position and momentum when it is in thermal equilibrium. This is not an issue classically, as we have many proofs outlining when systems are ergodic, meaning that the infinite time average is equal to the phase space average.

For quantum systems, however, such results are pretty much nonexistant. One of the first issues one has to deal with is that in a closed quantum system the evolution operator is unitary, meaning that not only is the energy conserved but the distance between two input states is preserved throughout the dynamics. A system $epsilon$ far away from thermal equilibrium at the beginning will remain $epsilon$ far away throughout the entire duration of the evolution. This is not what is observed in practice, as a quantum system placed in a dilution refrigerator will eventually cool down to the temperature of the fridge.

To resolve this apparent discrepancy physicists have settled on three main approaches: the Eigenstate Thermalization Hypothesis (ETH), Linbladian based evolutions, and the Repeated Interactions (RI) framework. ETH studies thermalization without giving up the notion of a closed quantum system. In this framework the system is taken to be large, as in thermodynamically large, and thermalization is observed whenever one considers the state of a single particle or a single _local_ observable. In chaotic systems these observables can be rigorously shown to appear as if they came from the thermal average $e^(-beta H_"local") / cal(Z)_"local"$, but showing that similar techniques work for non-chaotic systems is a major open question. In fact, the existance of many-body localization serves as a counterexample to a universal ETH, but the existance of these phases of matter in the real world, along with the validity of ETH in general, are heavily debated.

The second approach based on Linbladian schemes abandons the notion of a closed quantum system and instead studies only the effect of the environment on the system. This is a valid idea, as we know from basic quantum information theory that _any_ channel on a system can be mimicked with only a quadratically larger Hilbert space. This means that even the effects of an infinitely large environment can in principle be simulated on finite sized quantum devices. This is the approach taken by @davies1974markovian. Some of the downsides to this approach is that the specific models used for the environment typically make strong assumptions, such as weak coupling, Markovianity, or an infinite number of degrees of freedom, which can break down in specific scenarios.

The last approach is the Repeated Interactions framework which is somewhat of a hybrid of the two previous formalisms. The RI prescription does not abandon the notion of closed quantum systems but instead opts to directly simulate the environment via very small added degrees of freedom. The most straightforward example is to imagine a single photon $gamma$ drawn from some black-body thermal spectrum. This single photon then interacts with the system of interest for a brief period, only to then fly off and never interact with the system again. We then can take a new photon from this thermal background and interact with the system again, but with a refreshed photon. The idea is that by repeating this over and over, if the interaction is properly chosen, the many photons can collectively thermalize the system. The main drawback with this approach is that not only does a model for the environment need to be chosen, but physically realistic interactions need to be inserted in order for the system to converge to the thermal state. Existing RI results have shown that single spin environments are sufficient to thermalize single spin systems, or at the largest a three level system.

In this chapter we explore how to lift these restrictions on the RI framework through a randomized interaction approach. The resulting dynamics are not unitary but instead mixed unitary. We show how to analyze this channel in a weak-coupling regime, which is effectively a Taylor Series with respect to a coupling constant $alpha$. This expansion then reveals a Markov chain underlying the resulting channel on the system. By using basic tools from Markov Chain analysis we can then compute the fixed points of the map and bound the number of interactions needed to converge to the fixed point. By showing that the thermal states are the unique fixed points we can then prove thermalization. The beauty of our technique is that it works for any non-degenerate system with or without knowledge of the eigenvalues of the system, although without eigenvalue knowledeg we are only able to show that the thermal state is an approximate fixed point for finite $beta$ but exactly fixed in the limit $beta -> oo$.

The rest of this chapter is organized as follows. In @sec_tsp_intro we briefly discuss related works in quantum algorithms and provide a summary of the main technical results. In @sec_tsp_weak_coupling we develop the weak-coupling expansion and provide necessary analysis about the underlying Markov chain. In @sec_tsp_oscillator we then take these results and show how to use them to prepare thermal states for single qubit systems and for truncated harmonic oscillators. This section can be viewed as a warmup to the more general results contained in @sec_tsp_generic_sys, but we present the results separately as they utilize slightly different techniques and are more readily comparable to existing approaches. In @sec_tsp_generic_sys we study generic systems from two perspectives, one in which no eigenvalue knowledge is present and the other in which eigenvalues are known. Finally we include a discussion on interpretations of these results and possible extensions in @sec_tsp_discussion.

== Related Work and Main Results <sec_tsp_intro>



The simulation of quantum systems and materials is among the most promising applications for exponential advantages of digital quantum computers over classical computers @aspuru2005simulated @reiher2017elucidating @tensorHypercontraction. A critical step in quantum simulation algorithms, as well as other quantum algorithms such as Semi-Definite Program (SDP) solvers @brandao2019sdp and Hamiltonian learning routines @anshu_sample-efficient_2021, is the preparation of good input states, which are typically thermal states $frac(e^(-beta H), tr(e^(-beta H)))$. Thermal states at low temperatures (high $beta$) have large overlap with the ground states of the system, indicating that preparing thermal states is just as difficult as the QMA-Hard $k$-local ground state preparation problem @kempe2005complexitylocalhamiltonianproblem.

Many classical algorithms have been developed to estimate the measurement outcomes of quantum experiments with the workhorse behind many of these algorithms typically being some kind of Metropolis-Hastings algorithm @metropolis1953equation to implement a Markov Chain Monte Carlo (MCMC) program. The Metropolis-Hastings algorithm solves the problem of sampling from arbitrary probability distributions and can be used to estimate partition functions, a \#P-Hard problem @roth1996hardness. Despite the difficulty of the problem it solves and minimal theoretic guarantees on the runtime, the Metropolis-Hastings algorithm has worked resoundingly well in practice. This is in part due to its elegant simplicity and ease of implementation. The algorithm does tend to breakdown in a few important areas though, namely quantum systems with "sign problems" (@signProblemOG, @troyer2005sign), very high dimension systems @beskos2010optimaltuninghybridmontecarlo, and Hamiltonians with many deep local minima @betancourt2018conceptualintroductionhamiltonianmonte. The sign problem in particular serves as an important impetus for developing quantum computers, which naturally do not have to deal with it. However attempts to naively port the classical Metropolis-Hastings algorithm to quantum computers has been rather difficult due to inherent difficulties with quantum information, such as no-cloning. Initial attempts @temme2011 are rather cumbersome and attempts to deal with the filtering and rejection stages relying on "quantum unwinding" techniques developed by Marriott and Watrous @marriott2005quantum. These complications make the resulting quantum algorithms difficult to analyze.

In recent years new approaches have been developed @chen2023quantumthermalstatepreparation, @gilyen2024quantumgeneralizationsglaubermetropolis, @motlagh2024ground, @motta2019 @ding2024single, many of which are based on the simulation of Linblad operators from open quantum systems @davies1974markovian. These algorithms have seen a marked improvement in recent years, ranging from ground state preparation routines with single-ancilla overhead @ding2024single to the first constructions that satisfy a discrete-time detailed balance condition @gilyen2024quantumgeneralizationsglaubermetropolis. The correctness of many of these algorithms, such as @ding2024efficientquantumgibbssamplers, is based on satisfying the Kubo-Martin-Schwinger (KMS) condition (@kms2, @kms1) which guarantees that the thermal state is a fixed point of the dynamics. The literature on this class of algorithms is already significant and continues to grow, so we point the reader to @gilyen2024quantumgeneralizationsglaubermetropolis, @dalzell2023quantumalgorithmssurveyapplications, @chen2023quantumthermalstatepreparation and @rouze2024efficientthermalizationuniversalquantum for their thorough literature reviews.

One of the main drawbacks to the above approaches is the sheer complexity of the resulting algorithms. These algorithms tend to rely on coherently weighted sums of Heisenberg evolved jump operators and the construction of circuits to simulate the resulting Linbladians is nontrivial, as mentioned in Section 1.2 of @gilyen2024quantumgeneralizationsglaubermetropolis. Further, these algorithms tend to require logarithmically more ancilla qubits to allow for the addition of jump operators whereas our routine explicitly utilizes only a single ancilla qubit. Turning to ground states specifically, there exists single ancilla algorithms @ding2024single but we remark that our channel is the first general purpose thermal state preparation routine for finite $beta$ that utilizes only one qubit explicitly. Further, our routine avoids the complication of simulating weighted Linbladians and has incredibly simple circuits only relying on time independent Hamiltonian simulation and Haar 2-designs.

=== Main Results
The remainder of the paper is split into three main parts. Section @sec_tsp_weak_coupling contains a derivation of the weak-coupling expansion in Lemma [??] and outlines the underlying Markov chain behavior in Section [??]. This weak-coupling expansion is presented in as much generality as possible as it may be of use in other applications beyond our thermalization procedure, such as thermometry or spectroscopy. Section @sec_tsp_oscillator has two theorems concerning single qubit systems and harmonic oscillators, Theorems [??] and [??] respectively, as well as numerics exploring the $beta$ and $epsilon$ dependence of the channel. Section @sec_tsp_generic_sys contains our most general results in Theorems [??] and [??] in which we show that the thermal state is an approximate fixed point for arbitrary Hamiltonians, bound the runtime in terms of a Markovian spectral gap, and finally compute this spectral gap for the ground state limit. The main difference between these two theorems is that one requires only an uppper bound on the spectral norm $norm(H_S)$ while the other takes advantage of eigenvalue knowledge.

One of the key aspects of our thermalization procedure is that the analysis is dependent on the ability to tune the environment gap $gamma$ to match the system energy differences. One of our main results in Theorem [??] shows that even if the user cannot tune $gamma$ at all and is reduced to uniform guessing within an interval containing all the differences $Delta_S (i,j)$, then thermalization can still occur. We show that the thermal state at finite $beta$ is an approximate fixed state, with the error going to 0 as the coupling constant $alpha -> 0$. This zero coupling limit can be taken with the opposite limit $t -> infinity$ to yield a nonzero simulation time for the random interaction $G$. Further, we show that the ground state is exactly the fixed point in the $beta -> infinity$ limit. In this limit we are also able to bound the total simulation time required as $L dot t in tilde(O)(frac(dim_S^(16) norm(H_S)^7, delta_("min")^8 epsilon^6))$, where $delta_("min")$ represents a "resolution" type distance and is the smallest difference between two distinct eigenvalue differences $|Delta_S (i,j) - Delta_S (k,l)|$. When preparing finite $beta$ thermal states we pick up an extra factor of $frac(1, tilde(lambda)_star(beta)^7)$ related to the spectral gap of the transition matrix.

== Weak Coupling Expansion <sec_tsp_weak_coupling>

=== Preliminaries and Notation <sec_tsp_prelims>
We will be working with a bipartite Hilbert space consisting of a system space $hilb_S$ with dynamics governed by the Hamiltonian $H_S$ and an environment space $hilb_E$ with Hamiltonian $H_E$. The total space is $hilb = hilb_S tp hilb_E$ with Hamiltonian $H = H_S tp identity_E + identity_E tp H_E = H_S + H_E$. We will assume without loss of generality that our spaces are encoded in qubits so that $hilb_S = bb(C)^(2^n)$ and $hilb_E = bb(C)^(2^m)$. We use $dim_S$ to refer to the dimension of the system's Hilbert space ($2^n$), $dim_E$ the environment, and $dim$ the total Hilbert space. As for the basis we will use for our spaces, we will work directly in the eigenbasis of each Hamiltonian. Besides simplifying our calculations, we can do so because the interaction term we will introduce later is unitarily invariant. We denote these basis in a 1-indexed fashion as

$
    H_(S) = sum_(i = 1)^(2^n) lambda_S (i) ketbra(i, i) ,#h(1fr) H_(E) = sum_(j=1)^(2^m) lambda_E (j) ketbra(j, j) ,#h(1fr) H = sum_(i=1)^(2^n) sum_(j=1)^(2^m) lambda(i,j) ketbra(i\,j, i\,j),
$

where $lambda(i,j) = lambda_S(i) + lambda_E(j)$ and we will sort the eigenvalues in nondecreasing order such that $i > j => lambda_S(i) >= lambda_S(j)$. We note that the ground state in our 1-indexed notation is therefore $ketbra(1,1)$. We also make use of the following notation for the energy differences of the system-environment Hamiltonian and just the system

$ Delta(i,j|k,l) := lambda(i,j) - lambda(k,l), quad Delta_S(i,i') = lambda_S(i) - lambda_S(i'), $ <eq_delta_def>

and because our eigenvalues are sorted $i > j => Delta_S(i,j) >= 0$. We will need a few other notations for eigenvalue differences. First we denote the degeneracy of an eigenvalue $lambda(i,j)$ using $eta(i,j)$ and the number of times a system eigenvalue _difference_ is present as $eta_Delta (i,j)$. For example, in a truncated harmonic oscillator with 4 energy levels the lowest gap $Delta$ is present 3 times, so $eta_Delta (1, 2) = 3$. The second is that we will need to eventually analyze interferences between eigenvalue differences of the system, so we define

$ delta_(min) := min_(Delta_S(i,j) != Delta_S(k,l)) lr(| Delta_S(i,j) - Delta_S(k, l) |). $ <eq_delta_min_def>

Note that nothing in this definition prevents one of the summands, say $Delta_S (k,l)$, from being 0. This implies that $delta_(min) <= Delta_S (i,j)$ for all $i$ and $j$.

Currently our dynamics involved a system separated from the environment, so we need to fix this by adding an interaction term $G : hilb_S := hilb_E -> hilb_S tp hilb_E$. We will choose $G$ randomly via the eigendecomposition

$
    G = U_("haar") D U_("haar")^dagger, U_("haar") tilde "Haar"(hilb_S tp hilb_E) text(" and ") D_(i i) tilde cal(N)(0,1),
$ <eq_interaction_def>

where the eigenvectors are Haar distributed and the eigenvalues I.I.D. normal Gaussian variables. We then add this random interaction term to our system-environment dynamics with a coupling constant $alpha$, yielding a total dynamics governed by $H_S + H_E + alpha G$. We define the following rescaled coupling constant

$ tilde(alpha) := frac(alpha t, sqrt(dim + 1)), $ <eq_a_tilde_def>

where the $dim$ is the total Hilbert space $hilb$ dimension. The rescaling with respect to $dim$ is to capture the factors of $1/(dim + 1)$ in the transition amplitudes that appear later and leads to much more compact expressions.
This gives a decomposition of expectation values over $G$ into two parts

$ EE_G f(G) = EE_("haar") EE_(D) f(G), $

where the two expectations on the right commute with each other $bb(E)_("haar") bb(E)_(D) = bb(E)_(D) bb(E)_("haar")$.

We will use this interaction term to couple our system to an environment prepared in the thermal state $rho_E(beta) = e^(-beta H_E) /partfun_E(beta)$, where $partfun_E(beta) = tr(e^(-beta H_E))$, and then trace out the environment. This gives the definition of our thermalizing channel $Phi : cal(L)(hilb_S) -> cal(L)(hilb_S)$ as

$
    Phi(rho \; alpha, beta, t) := tr_(hilb_E) bb(E)_(G) lr([ e^(+i(H + alpha G)t) rho tp rho_E(beta) e^(-i(H + alpha G) t)]) .
$ <eq:PhiDef>

We will typically drop the implicit parameters of $alpha, beta$ and $t$. Our goal is to show how this channel can be used to prepare the system in the thermal state $rho(beta) = frac(e^(-beta H_S), partfun(beta))$. It will be useful to introduce a fixed-interaction channel $Phi_G : cal(L)(hilb_S tp hilb_E) -> cal(L)(hilb_S tp hilb_E)$ over the total Hilbert space $hilb$ as

$
    Phi_G (rho tp rho_E; alpha, t) := e^(+i(H + alpha t)) rho tp rho_E e^(- i(H + alpha G)t),
$ <eq_phi_g_definition>

giving us $Phi (rho\; alpha, beta, t) = tr_(hilb_E) bb(E)_G Phi_G (rho tp rho_E(beta); alpha, t)$. Another alternative notation for $Phi$ that we will use is whenever $hilb_E$ is a single qubit with energy gap $gamma$ we will use $Phi_gamma$ to draw attention to this specific energy gap. We will also make frequent use of indicator functions, denoted $bold(I)[P]$, which is 1 if the predicate $P$ is true and 0 if $P$ is false.

=== First and Second Order Expansion <sec_tsp_expansion_series>

In order to understand our thermalizing channel $Phi$ we will compute a Taylor Series for the output of the channel with respect to the coupling constant $alpha$. We will perform the $alpha$ expansion about $alpha = 0$ and we will use the mean value form of the remainder, in which we are guaranteed a special value $alpha_(star) in (0, infinity)$ such that the final derivative evaluated at $alpha_(star)$ is the exact amount needed. We use a second-order expansion and will need to explicitly compute terms up to order $alpha^2$, which will give the following expansion
$
    Phi (rho; alpha) = Phi (rho; 0) + alpha frac(partial, partial alpha) Phi (rho; alpha) |_(alpha = 0) + frac(alpha^2, 2) frac(partial^2, partial alpha^2) Phi (rho; alpha) |_(alpha = 0) + R_(Phi)(rho; alpha_(star)).
$ <eq_tsp_phi_taylor_series>
We use
$
    cal(T)(rho) := frac(alpha^2, 2) frac(partial^2, partial alpha^2) Phi (rho; alpha) |_(alpha = 0) = frac(alpha^2, 2) tr_(hilb_E) bb(E)_(G) lr([frac(partial^2, partial alpha^2) Phi_G(rho; alpha) bar.v _(alpha = 0)])
$ <eq:transition_def>

to denote the transition terms, as it will be revealed that the first two terms do not cause transitions in the system state, and $R_(Phi)$ to denote the remainder.
Further we will often leave the dependence on the $alpha$ parameter implicit and only include it when necessary.

We start off with the $O(alpha^0)$ term, which can be trivially computed as
$
    Phi (rho; 0) = tr_(hilb_E) integral e^(i(H + alpha G) t) rho tp rho_E(beta) e^(-i (H + alpha G) t) d G |_(alpha = 0) = e^(i H t) rho e^(-i H t).
$

We then see that if $[ rho, H] = 0$ then $Phi (rho; 0) = id (rho)$, and as we restrict ourselves to such input states we will use this throughout the remainder of the paper. The next order correction is the $O(alpha^1)$ term.
#theorem([First Order $Phi$])[
    Let $Phi$ be the thermalizing quantum channel given by @eq:PhiDef and $G$ the randomly chosen interaction term as given by @eq_interaction_def. The $O(alpha)$ term in the weak-coupling expansion in @eq_tsp_phi_taylor_series vanishes
    $ frac(partial, partial alpha) Phi (rho; alpha) |_(alpha = 0) = 0. $
] <thm_tsp_first_order_phi>
#proof()[
    We start by moving the $alpha$ derivative through the linear operations of partial tracing and integrals so that it can act on the fixed interaction map $Phi_G$
    $
        frac(partial, partial alpha) Phi (rho) |_(alpha = 0) &= frac(partial, partial alpha) tr_(cal(H)_E)(integral Phi_G (rho) d G) |_(alpha = 0) \
        &= tr_(cal(H)_E)(integral frac(partial, partial alpha) Phi_G(rho) d G |_(alpha = 0)) .
    $
    Now we use the expression for $Phi_G$ in Eq. @eq_phi_g_definition to compute the derivatives,
    $
        frac(partial, partial alpha) Phi_G (rho) =& (frac(partial, partial alpha) e^(+ i (H + alpha G)t)) rho tp rho_E e^(-i (H + alpha G) t) + e^(+i (H + alpha G)t) rho tp rho_E (frac(partial, partial alpha) e^(- i (H + alpha G)t)) \
        =& (integral_(0)^(1) e^(i s (H+alpha G)t) (i t G) e^(i (1-s) (H+alpha G)t) d s) rho tp rho_E e^(-i(H+alpha G)t) \
        &+ e^(i(H+alpha G)t) rho tp rho_E (integral_(0)^1 e^(-i s (H+alpha G) t) (- i t G) e^(-i (1-s) (H+alpha G)t) d s).
    $ <eq_first_order_alpha_derivative>
    We can set $alpha = 0$ in the above and introduce the expectation over $G$ that will be required
    $
        bb(E)_G[ frac(partial, partial alpha) Phi_G(rho) |_(alpha = 0)] &= i t bb(E)_G integral_0^1 e^(i s H t) G e^(-i s H t) d s e^(i H t) rho tp rho_E e^(-i H t) \
        &- i t e^(+i H t) rho tp rho_E bb(E)_G integral_0^1 e^(-i s H t) G e^(-i(1-s) H t) d s \
        &= i t integral_0^1 e^(i s H t) bb(E)_G [G] e^(-i s H t) d s e^(i H t) rho tp rho_E e^(-i H t) \
        &- i t e^(+i H t) rho tp rho_E integral_0^1 e^(-i s H t) bb(E)_G [G] e^(-i(1-s) H t) d s.
    $
    Since our eigenvalues, $D_(i i)$, are mean zero ($EE_D D = 0$) we can compute $bb(E)_G [G]$ and arrive at the lemma statement
    $
        EE_G [G] = EE_"haar" EE_D [U_"haar" D U_"haar"^dagger] = EE_"haar" [U_"haar" EE_D [D] U_"haar"^dagger] = 0.
    $
]

#h(5mm) Now we move on to the $O(alpha^2)$ term in the weak-coupling expansion of $Phi$. We first will compute the combined system-environment output of a generic system-environment basis state and we note that this result holds for an arbitrary dimension environment. We will use this to draw two results: the first being for a single qubit environment the transition amplitudes of just the system can be split into on-resonance and off-resonance terms based on the tuning of the environment qubit Hamiltonian. The second result is that coherences are not introduced to the state at this order of $Phi$, meaning if an input density matrix $rho$ is diagonal then $(id + cal(T))(rho)$ will also be diagonal. This will be crucial for our later understanding of the channel as a Markov chain.
#lemma()[
    Given a system Hamiltonian $H_S$, an environment Hamiltonian $H_E$, a simulation time $t$, and coupling coefficient $alpha$, let $Phi_G$ denote the time evolution channel under a fixed interaction term $G$ as given in @eq_phi_g_definition, let $chi$ denote the following coherence prefactor
    $
        chi(i,j) := sum_(a,b: Delta (i,j,|a,b) != 0) (1 - i Delta (i,j|a,b)t - e^(-i Delta (i,j|a,b) t)) / (Delta (i,j|a,b)^2),
    $
    and let $eta(i,j)$ denote the degeneracy of the $(i,j)^"th"$ eigenvalue of $H = H_S + H_E$. Then the $O(alpha^2)$ term of $Phi_G$ in a weak-coupling expansion is given by
    $
        &alpha^2 / 2 EE_G [diff^2 / (diff alpha^2) Phi_G (ketbra(i \,j, k \,l))|_(alpha = 0) ] \
        =& - (alpha^2 e^(i Delta (i,j|k,l) t)) / (dim + 1) (chi(i,j) + chi(k,l)^* + t^2 / 2 (eta(i,j) + eta(k,l)) )ketbra(i\,j, k\,l) \
        &+ angle.l i,j|k,l angle.r (alpha^2 t^2) / (dim + 1) sum_(a,b) sinc^2 ( Delta(i,j|a,b) t / 2) ketbra(a\,b, a\,b)
    $ <eq_el_gigante>
    For $ket(i \, j) = ket(k \, l)$ the above expression simplifies to
    $
        &alpha^2 / 2 EE_G [diff^2 / (diff alpha^2) Phi_G (ketbra(i \,j, i \,j))|_(alpha = 0) ] \
        &= tilde(alpha)^2 sum_((a,b) != (i,j)) sinc^2 (Delta (i,j|a,b) t / 2) (ketbra(a\, b, a\, b) - ketbra(i\, j, i \,j))
    $ <eq_el_gigante_dos>
    which also demonstrates that $tr cal(T)(rho) = 0$ for $rho$ such that $[rho, H_S] = 0$.
] <lem_tsp_transitions>
The proof of this lemma uses similar techniques to the proof of @thm_tsp_first_order_phi but is significantly more technical and can be found in @sec_appendix_haar.


Next we will compute the effects of the channel on just the system alone which involves computing the partial trace $tr_(hilb_E)$. We can either do this for a generic environment, which results in summations over $hilb_E$ floating around, or specialize to a specific choice of $hilb_E$ and compute the summation. For the remainder of this paper we will choose the latter option with a single qubit environment $hilb_E = CC^2$ and denote the Hamiltonian $H_E = mat(0, 0; 0, gamma)$. Our environment input states then become
$
    rho_E (beta) = (e^(-beta H_E)) / (cal(Z)_E (beta)) = 1 / (1 + e^(-beta gamma)) ketbra(0, 0) + e^(-beta gamma) / (1 + e^(-beta gamma)) ketbra(1,1) := q(0) ketbra(0,0) + q(1) ketbra(1,1) ,
$ <eq_env_state_def>
where we will use the environment qubit probabilities $q(0)$ and $q(1)$ in calculations for brevity. It will turn out that the value chosen for $gamma$ is highly critical to the convergence of our algorithm, tuning it to match eigenvalue _differences_ of the system $H_S$ will allow us to analyze the convergence of the algorithm. As we can see in @eq_el_gigante there will be a lot of $sinc$ functions used, we will characterize a $sinc$ function as being on-resonance or off-resonance if the inputs are sufficiently close to zero (the max for sinc). As for how close "sufficiently close" actually is will depend on various parameters, such as $t, alpha, epsilon$, and the spectral properties of $H_S$.
#theorem([Second-Order Expansion $cal(T)$])[
    Let $cal(T)$ denote the second-order correction for a weak coupling expansion for a thermalizing channel $Phi$ with a single qubit environment. The following properties hold for the second order term.
    + The transition element from $ketbra(i,i)$ to $ketbra(j,j)$, for $i != j$, is given by $
    bra(j) cal(T) (ketbra(i,i)) ket(j) = tilde(alpha)^2 (&sinc^2(Delta_S (i,j) t/2) + 1/(1 + e^(-beta gamma)) sinc^2 ((Delta_S (i,j) - gamma)t/2) \
    & + e^(-beta gamma)/(1 + e^(-beta gamma)) sinc^2((Delta_S (i,j) + gamma) t/2)).
  $ <eq_transition_terms_total>
    + For same-state transitions $ketbra(i,i)$ to $ketbra(i,i)$ we have $
    bra(i) cal(T)(ketbra(i,i)) ket(i) = - sum_(j != i) bra(j) cal(T) (ketbra(i,i)) ket(j),
  $ which follows from $tr cal(T)(rho) = 0$ as shown in @lem_tsp_transitions.
    + There are no coherences, or off-diagonal density matrix elements, introduced in the system up to $O(alpha^2)$, or mathematically $ j != k ==> bra(j) cal(T) (ketbra(i,i)) ket(k) = 0. $
]<thm_tsp_second_order_expansion>
Before we prove this result we will introduce the concept of on- and off-resonant transitions which we give below.
#definition([On and Off Resonant Transitions])[
    The transition elements in @eq_transition_terms_total can be divided into on-resonance and off-resonance transitions based on the arguments to the $sinc$ function. We define the on-resonance transitions as
    $
        bra(j) cal(T)_"on" (ketbra(i,i)) ket(j) := & tilde(alpha)^2 q0 II[ |Delta_S (i,j) - gamma| <= delta_min] sinc^2 ((Delta_S (i,j) - gamma ) t / 2) \
        + & tilde(alpha)^2 q1 II[ |Delta_S (i,j) + gamma| <= delta_min] sinc^2 ((Delta_S (i,j) + gamma) t / 2)
    $ <eq_on_resonance>
    and the off-resonance terms as
    $
        bra(j) cal(T)_"off" (ketbra(i,i)) ket(j) := & tilde(alpha)^2 q0 II[ |Delta_S (i,j) - gamma| > delta_min] sinc^2 ((Delta_S (i,j) - gamma ) t / 2) \
        + & tilde(alpha)^2 q1 II[ |Delta_S (i,j) + gamma| > delta_min] sinc^2 ((Delta_S (i,j) + gamma) t / 2) \
        + & tilde(alpha)^2 sinc^2 (Delta_S (i,j) t / 2).
    $ <eq_off_resonance>
    For the same-state transitions $ketbra(i,i)$ to $ketbra(i,i)$ the on- and off-resonance transitions are equal to
    $
        bra(i) cal(T)_"on" (ketbra(i,i)) ket(i) &= - sum_(j != i) bra(j) cal(T)_"on" (ketbra(i,i)) ket(j) \
        " and " bra(i) cal(T)_"off" (ketbra(i,i)) ket(i) &= - sum_(j != i) bra(j) cal(T)_"off" (ketbra(i,i)) ket(j).
    $ <eq_same_state_transition_resonances>
] <def_transition>
We will now use these definitions to prove @thm_tsp_second_order_expansion.
#proof([of @thm_tsp_second_order_expansion])[
    The bulk of this proof will be based on straightforward reductions from @eq_el_gigante. To start we will first show that no off-diagonal elements are introduced to the density matrix. By taking the $(j,k)$ matrix element of the output from @eq_el_gigante we see
    $
        bra(j) cal(T)(ketbra(i,i))ket(k) &= sum_(l, m) e^(-beta lambda_E (m)) / (1 + e^(-beta lambda_E (m))) bra(j\, l) alpha^2 / 2 EE_G [ partial^2 / (diff alpha^2) Phi_G (ketbra(i\, m, i\, m)) |_(alpha = 0) ] ket(k\, l) \
        &= - sum_(l,m) q(m) tilde(alpha)^2 (chi (i,m) + chi (i,m)^* + t^2 eta (i,m)) braket(j\, l, i\, m) braket(i\, m, k\, l) \
        & + sum_(l, m) q(m) sum_(a,b) tilde(alpha)^2 sinc^2 (Delta (i,m|a,b) t / 2 ) braket(j\,l, a\,b) braket(a\, b, k\, l) \
        &= 0,
    $
    where we introduce $q(m)$ for $m=0,1$ to be a placeholder for the prefactors in @eq_el_gigante and the last equality is due to the fact that $j != k$ implies that $braket(j\, l, i\,m)$ and $braket(i\,m, k\, l)$ cannot both be nonzero and likewise for $braket(j\, l, a\,b)$ and $braket(a\,b, k\,l)$.

    Since we have shown that coherences are not introduced to our system we can focus on the transitions from diagonal entries to diagonal entries in $rho$. We make heavy use of Eq. \eqref{eq:el_gigante_dos} which tells us that for $i != k$ the system-environment transition amplitude is
    $
        alpha^2 / 2 bra(k\, l) EE_G [ diff^2 / (diff alpha^2) Phi_G (ketbra(i\, j, i\,j)) |_(alpha = 0) ] ket(k\, l) = tilde(alpha)^2 sinc^2 ( Delta (i,j | k, l) t / 2) .
    $
    Now because all the operations present in the above expression are linear we can compute this map for the initial environment state $rho_E (beta)$ straightforwardly. Taking the output of this linear combination and computing the trace over the environment then gives us the expression for $cal(T)$ using the assumption that the environment is a single qubit we find using the definition of $gamma$ and $Delta_S$ in @eq_delta_def
    $
        bra(j) cal(T)(ket(i)bra(i)) ket(j) &= sum_(k, l) q(k) frac(alpha^2, 2)bra(j\, l) bb(E)_G [frac(partial^2, partial alpha^2) Phi_G(ket(i\, k)bra(i\,k)) |_(alpha = 0)] ket(j\, l) \
        &= tilde(alpha)^2 sum_(k, l) q(k) sinc^2 (frac(Delta(i, k | j , l) t, 2)) \
        &= tilde(alpha)^2 (q(0) sinc^2 (frac(Delta(i, 0 | j , 0) t, 2)) + q(0) sinc^2 (frac(Delta(i, 0 | j , 1) t, 2))) \
        & quad + tilde(alpha)^2 (q(1) sinc^2 (frac(Delta(i, 1 | j , 0) t, 2)) + q(1) sinc^2 (frac(Delta(i, 1 | j , 1) t, 2))) \
        &= tilde(alpha)^2 (q(0) sinc^2 (frac(Delta_S (i,j) t, 2)) + q(0) sinc^2 (frac((Delta_S (i,j) - gamma) t, 2))) \
        & quad + tilde(alpha)^2 (q(1) sinc^2 (frac((Delta_S (i, j) + gamma) t, 2)) + q(1) sinc^2 (frac(Delta_S (i,j) t, 2))),
    $
    where we see that combining the terms with $sinc^2 (Delta_S(i,j) t / 2)$, as $q(0) + q(1) = 1$, we immediately get @eq_transition_terms_total.

    To classify these terms as on-resonance or off-resonance we will focus on the argument to the sinc function, which is of the form $Delta_S(i,j) t/ 2$ or $(Delta_S(i,j) plus.minus gamma) t/ 2$. The idea is that we will take $t$ large enough so that only the energy differences that are less than $delta_"min"$, as defined in @eq_delta_min_def, will be non-negligible. Clearly the term $tilde(alpha)^2 sinc^2 ( frac(Delta_S (i,j)t, 2) )$ will always be off-resonance, as $delta_"min" <= Delta_S (i,j)$.

    Now we have three terms to classify as either on-resonance or off-resonance, we refer to each term by their argument to the $sinc$ function. The first we can categorically declare as being off-resonance is the $Delta_S(i,j)$ term. By [??] we know $sinc^2(Delta_S(i,j) t/ 2) <= 4 / (delta_"min"^2 t^2)$, which we will make arbitrarily small in later sections. The other two can only be classified as on or off resonance depending if $Delta_S(i,j)$ is positive or negative. If $i > j$ then we know that $Delta_S(i,j) >= 0$ and therefore $sinc^2((Delta_S(i,j) - gamma)t/2)$ term can be close to 1 if $gamma approx Delta_S(i,j)$, which also shows the $Delta_S(i,j) + gamma$ term is off-resonance for all $gamma$. We say that the $Delta_S(i,j) - gamma$ term in this scenario is on-resonance if $|Delta_S(i,j) - gamma| <= delta_"min"$. This classification is best described symbolically as
    $
        i > j "and" |Delta(i,j) - gamma| <= delta_("min") ==>
        bra(j) cal(T)_"on" (ketbra(i,i)) ket(j) =
        tilde(alpha)^2 q(0) "sinc"^2 ( (Delta_S (i,j) - gamma)t / 2 ).
    $ <tmp_tsp_1>
    The $q(0)$ prefactor indicates that the ancilla started in it's low energy state and since $sinc^2$ is symmetric we can write the argument as $gamma Delta_S (i,j)$ which can be remembered as the ancilla gaining $gamma$ amount of energy and the system losing $Delta_S (i,j)$. In this scenario the $Delta_S (i,j) + gamma$ term is therefore put in the off-resonance map
    $
        i > j "and" |Delta_S (i,j) - gamma| <= delta_min \
        ==>
        bra(j) cal(T)_"off" (ketbra(i,i)) ket(j) =
        tilde(alpha)^2 ( sinc^2 ( Delta_S (i,j) t / 2 ) +
            q(1) sinc^2 ( (Delta_S (i,j) + gamma) t / 2 ) ).
    $

    Now for $i < j$ we find that the on-resonance term is
    $
        i < j "and" |Delta_S (i,j) + gamma| <= delta_("min") ==>
        bra(j) cal(T)_"on" (ketbra(i, i)) ket(j) =
        tilde(alpha)^2 q(1) sinc^2 ( (Delta_S (i,j) + gamma)t / 2 ).
    $ <tmp_tsp_2>
    Similarly to before the $q(1)$ prefactor tells us the ancilla starts in the excited state. This matches with the energy argument by noting that $Delta_S (i,j) <= 0$ and that the argument to $"sinc"$ is symmetric, which allows us to write it as $|Delta_S (i,j)| - gamma$; indicating that the system gains energy $|Delta_S (i,j)|$ and the ancilla energy _drops_ by $-gamma$ (therefore increases by $gamma$). In this scenario the $Delta_S (i,j) - gamma$ term is off-resonance and we have
    $
        i < j "and" |Delta_S (i,j) + gamma| <= delta_("min") \
        ==>
        bra(j) cal(T)_"off" (ketbra(i, i)) ket(j) =
        tilde(alpha)^2 ( sinc^2 ( Delta_S (i,j) t / 2 ) +
            q(0) sinc^2 ( (Delta_S (i,j) - gamma) t / 2 ) )
    $

    Now to compute the $i = j$ case, it is sufficient to utilize our results from the $i != j$ scenario. This is because our second order correction has zero trace $tr(cal(T)(rho)) = 0$ from , so we can define the on-resonance and off-resonance terms as the following

    $
        bra(i) cal(T)(ket(i)bra(i)) ket(i) &= - tilde(alpha)^2 sum_(k != i) bra(k) cal(T)(ketbra(i,i)) ket(k) \
        &= - tilde(alpha)^2 sum_(k != i) bra(k) (cal(T)_"on" (ketbra(i, i)) + cal(T)_"off" (ket(i)bra(i))) ket(k) \
        &=: bra(i) cal(T)_"on" (ketbra(i, i)) ket(i) + bra(i) cal(T)_"off" (ket(i)bra(i)) ket(i)
    $

    By plugging in @tmp_tsp_1 and @tmp_tsp_2 to the above we are done with the self-transition terms.
]

=== Markovian Dynamics and Error Terms <sec_tsp_markovian_dynamics>

Now that we have fully computed the significant contributors to the output of our channel $Phi$, we move on to characterize the behavior of the channel as a Markov chain with noise.
A Markov chain is a random process that involves a walker transitioning to vertices on a graph wherein the probability of transition does not depend on the history of the walker. Specifically, in this context we view the vertices in this graph as the eigenstates of the Hamiltonian. The repeated interaction model because of the lack of coherences in the weak coupling limit can be interpreted as a Markov process over these eigenstates with transitions probabilities given by the above analysis.

Specifically, the Markov chain is dictated by the $Phi (rho; 0)$ and $cal(T)_"on"$ terms in the weak-coupling expansion, for $[rho, H_S] = 0$ we showed that $Phi (rho; 0) = id (\ho)$, so from now on we will specifically only deal with such density matrices and characterize the zeroth order term as an identity map. As for the Markov chain, we will use normal font to denote matrices, such as $I$ for the identity matrix and $T$ for the transition term added on. We use $e_i$ to denote the basis vector associated with the quantum state $ketbra(i, i)$ and $p$ to denote the probability vector for $rho$ associated with its eigenvalues.


== Single Qubit and Truncated Harmonic Oscillator <sec_tsp_oscillator>

== Generic Systems <sec_tsp_generic_sys>

== Discussion <sec_tsp_discussion>

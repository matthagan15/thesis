#import "conf.typ": *

#set math.equation(number-align: bottom)

#let q0 = $1 / (1 + e^(-beta gamma))$
#let q1 = $e^(-beta gamma) / (1 + e^(-beta gamma))$

#import "@preview/ctheorems:1.1.3": *
#show: thmrules.with(qed-symbol: $square$)

#heading("Preparing Thermal Quantum States", level: 1, supplement: "Chapter") <ch:thermal_state_prep>
Thermal states of the form $e^(-beta H) / cal(Z)$ are ubiquitous in physics. These are the states that we believe physical systems take whenever they are cooled (or heated) to an inverse temperature of $beta$. They are typically used to estimate observables $O$ of interest in physically relevant states, with the thermal average typically denoted $angle.l O angle.r_beta = tr(O e^(-beta H) / cal(Z))$. These observables could be anything from dipole moments, magnetizations, or two body correlators. Classically these states correspond to the canonical ensemble, a distribution over phase space that tells us the probability of finding a particle at a particular position and momentum when it is in thermal equilibrium. This is not an issue classically, as we have many proofs outlining when systems are ergodic, meaning that the infinite time average is equal to the phase space average.

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

One of the key aspects of our thermalization procedure is that the analysis is dependent on the ability to tune the environment gap $gamma$ to match the system energy differences. One of our main results in Theorem [??] shows that even if the user cannot tune $gamma$ at all and is reduced to uniform guessing within an interval containing all the differences $Delta_S (i,j)$, then thermalization can still occur. We show that the thermal state at finite $beta$ is an approximate fixed state, with the error going to 0 as the coupling constant $alpha -> 0$. This zero coupling limit can be taken with the opposite limit $t -> infinity$ to yield a nonzero simulation time for the random interaction $G$. Further, we show that the ground state is exactly the fixed point in the $beta -> infinity$ limit. In this limit we are also able to bound the total simulation time required as $L dot t in tilde(O)(frac(dim_S^(16) norm(H_S)^7, delta_("min")^8 epsilon^6))$, where $delta_("min")$ represents a "resolution" type distance and is the smallest difference between two distinct eigenvalue differences $|Delta_S (i,j) - Delta_S (k,l)|$. When preparing finite $beta$ thermal states we pick up an extra factor of $frac(1, tilde(lambda)_star (beta)^7)$ related to the spectral gap of the transition matrix.

== Weak Coupling Expansion <sec_tsp_weak_coupling>

=== Preliminaries and Notation <sec_tsp_prelims>
We will be working with a bipartite Hilbert space consisting of a system space $hilb_S$ with dynamics governed by the Hamiltonian $H_S$ and an environment space $hilb_E$ with Hamiltonian $H_E$. The total space is $hilb = hilb_S tp hilb_E$ with Hamiltonian $H = H_S tp identity_E + identity_E tp H_E = H_S + H_E$. We will assume without loss of generality that our spaces are encoded in qubits so that $hilb_S = bb(C)^(2^n)$ and $hilb_E = bb(C)^(2^m)$. We use $dim_S$ to refer to the dimension of the system's Hilbert space ($2^n$), $dim_E$ the environment, and $dim$ the total Hilbert space. As for the basis we will use for our spaces, we will work directly in the eigenbasis of each Hamiltonian. Besides simplifying our calculations, we can do so because the interaction term we will introduce later is unitarily invariant. We denote these basis in a 1-indexed fashion as

$
    H_(S) = sum_(i = 1)^(2^n) lambda_S (i) ketbra(i, i) ,#h(1fr) H_(E) = sum_(j=1)^(2^m) lambda_E (j) ketbra(j, j) ,#h(1fr) H = sum_(i=1)^(2^n) sum_(j=1)^(2^m) lambda(i, j) ketbra(i\,j, i\,j),
$

where $lambda(i, j) = lambda_S (i) + lambda_E(j)$ and we will sort the eigenvalues in nondecreasing order such that $i > j => lambda_S (i) >= lambda_S (j)$. We note that the ground state in our 1-indexed notation is therefore $ketbra(1, 1)$. We also make use of the following notation for the energy differences of the system-environment Hamiltonian and just the system

$ Delta(i, j|k, l) := lambda(i, j) - lambda(k, l), quad Delta_S(i,i') = lambda_S (i) - lambda_S (i'), $ <eq_delta_def>

and because our eigenvalues are sorted $i > j => Delta_S(i,j) >= 0$. We will need a few other notations for eigenvalue differences. First we denote the degeneracy of an eigenvalue $lambda(i, j)$ using $eta(i, j)$ and the number of times a system eigenvalue _difference_ is present as $eta_Delta (i,j)$. For example, in a truncated harmonic oscillator with 4 energy levels the lowest gap $Delta$ is present 3 times, so $eta_Delta (1, 2) = 3$. The second is that we will need to eventually analyze interferences between eigenvalue differences of the system, so we define

$ delta_(min) := min_(Delta_S(i,j) != Delta_S(k,l)) lr(| Delta_S(i,j) - Delta_S(k, l) |). $ <eq_delta_min_def>

Note that nothing in this definition prevents one of the summands, say $Delta_S (k,l)$, from being 0. This implies that $delta_(min) <= Delta_S (i,j)$ for all $i$ and $j$.

Currently our dynamics involved a system separated from the environment, so we need to fix this by adding an interaction term $G : hilb_S := hilb_E -> hilb_S tp hilb_E$. We will choose $G$ randomly via the eigendecomposition

$
    G = U_("haar") D U_("haar")^dagger, U_("haar") tilde "Haar"(hilb_S tp hilb_E) text(" and ") D_(i i) tilde cal(N)(0,1),
$ <eq_interaction_def>

where the eigenvectors are Haar distributed and the eigenvalues I.I.D. normal Gaussian variables. We then add this random interaction term to our system-environment dynamics with a coupling constant $alpha$, yielding a total dynamics governed by $H_S + H_E + alpha G$. We define the following rescaled coupling constant

$ tilde(alpha) := frac(alpha t, sqrt(dim + 1)), $ <eq_a_tilde_def>

where the $dim$ is the total Hilbert space $hilb$ dimension. The rescaling with respect to $dim$ is to capture the factors of $1 / (dim + 1)$ in the transition amplitudes that appear later and leads to much more compact expressions.
This gives a decomposition of expectation values over $G$ into two parts

$ EE_G f(G) = EE_("haar") EE_(D) f(G), $

where the two expectations on the right commute with each other $bb(E)_("haar") bb(E)_(D) = bb(E)_(D) bb(E)_("haar")$.

We will use this interaction term to couple our system to an environment prepared in the thermal state $rho_E(beta) = e^(-beta H_E) / partfun_E(beta)$, where $partfun_E(beta) = tr(e^(-beta H_E))$, and then trace out the environment. This gives the definition of our thermalizing channel $Phi : cal(L)(hilb_S) -> cal(L)(hilb_S)$ as

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
    cal(T)(rho) := frac(alpha^2, 2) frac(partial^2, partial alpha^2) Phi (rho; alpha) |_(alpha = 0) = frac(alpha^2, 2) tr_(hilb_E) bb(E)_(G) lr([frac(partial^2, partial alpha^2) Phi_G(rho; alpha) bar.v_(alpha = 0)])
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
        chi(i, j) := sum_(a,b: Delta (i,j,|a,b) != 0) (1 - i Delta (i,j|a,b)t - e^(-i Delta (i,j|a,b) t)) / (Delta (i,j|a,b)^2),
    $
    and let $eta(i, j)$ denote the degeneracy of the $(i,j)^"th"$ eigenvalue of $H = H_S + H_E$. Then the $O(alpha^2)$ term of $Phi_G$ in a weak-coupling expansion is given by
    $
        &alpha^2 / 2 EE_G [diff^2 / (diff alpha^2) Phi_G (ketbra(i \,j, k \,l))|_(alpha = 0) ] \
        =& - (alpha^2 e^(i Delta (i,j|k,l) t)) / (dim + 1) (chi(i, j) + chi(k, l)^* + t^2 / 2 (eta(i, j) + eta(k, l)) )ketbra(i\,j, k\,l) \
        &+ angle.l i,j|k,l angle.r (alpha^2 t^2) / (dim + 1) sum_(a,b) sinc^2 ( Delta(i, j|a, b) t / 2) ketbra(a\,b, a\,b)
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
    rho_E (beta) = (e^(-beta H_E)) / (cal(Z)_E (beta)) = 1 / (1 + e^(-beta gamma)) ketbra(0, 0) + e^(-beta gamma) / (1 + e^(-beta gamma)) ketbra(1, 1) := q(0) ketbra(0, 0) + q(1) ketbra(1, 1) ,
$ <eq_env_state_def>
where we will use the environment qubit probabilities $q(0)$ and $q(1)$ in calculations for brevity. It will turn out that the value chosen for $gamma$ is highly critical to the convergence of our algorithm, tuning it to match eigenvalue _differences_ of the system $H_S$ will allow us to analyze the convergence of the algorithm. As we can see in @eq_el_gigante there will be a lot of $sinc$ functions used, we will characterize a $sinc$ function as being on-resonance or off-resonance if the inputs are sufficiently close to zero (the max for sinc). As for how close "sufficiently close" actually is will depend on various parameters, such as $t, alpha, epsilon$, and the spectral properties of $H_S$.
#theorem([Second-Order Expansion $cal(T)$])[
    Let $cal(T)$ denote the second-order correction for a weak coupling expansion for a thermalizing channel $Phi$ with a single qubit environment. The following properties hold for the second order term.
    + The transition element from $ketbra(i, i)$ to $ketbra(j, j)$, for $i != j$, is given by $ bra(j) cal(T) (ketbra(i, i)) ket(j) = tilde(alpha)^2 (&sinc^2(Delta_S (i,j) t / 2) + 1 / (1 + e^(-beta gamma)) sinc^2 ((Delta_S (i,j) - gamma)t / 2) \
            & + e^(-beta gamma) / (1 + e^(-beta gamma)) sinc^2((Delta_S (i,j) + gamma) t / 2)). $ <eq_transition_terms_total>
    + For same-state transitions $ketbra(i, i)$ to $ketbra(i, i)$ we have $ bra(i) cal(T)(ketbra(i, i)) ket(i) = - sum_(j != i) bra(j) cal(T) (ketbra(i, i)) ket(j), $ which follows from $tr cal(T)(rho) = 0$ as shown in @lem_tsp_transitions.
    + There are no coherences, or off-diagonal density matrix elements, introduced in the system up to $O(alpha^2)$, or mathematically $ j != k ==> bra(j) cal(T) (ketbra(i, i)) ket(k) = 0. $
] <thm_tsp_second_order_expansion>
Before we prove this result we will introduce the concept of on- and off-resonant transitions which we give below.
#definition([On and Off Resonant Transitions])[
    The transition elements in @eq_transition_terms_total can be divided into on-resonance and off-resonance transitions based on the arguments to the $sinc$ function. We define the on-resonance transitions as
    $
        bra(j) cal(T)_"on" (ketbra(i, i)) ket(j) := & tilde(alpha)^2 q0 II[ |Delta_S (i,j) - gamma| <= delta_min] sinc^2 ((Delta_S (i,j) - gamma ) t / 2) \
        + & tilde(alpha)^2 q1 II[ |Delta_S (i,j) + gamma| <= delta_min] sinc^2 ((Delta_S (i,j) + gamma) t / 2)
    $ <eq_on_resonance>
    and the off-resonance terms as
    $
        bra(j) cal(T)_"off" (ketbra(i, i)) ket(j) := & tilde(alpha)^2 q0 II[ |Delta_S (i,j) - gamma| > delta_min] sinc^2 ((Delta_S (i,j) - gamma ) t / 2) \
        + & tilde(alpha)^2 q1 II[ |Delta_S (i,j) + gamma| > delta_min] sinc^2 ((Delta_S (i,j) + gamma) t / 2) \
        + & tilde(alpha)^2 sinc^2 (Delta_S (i,j) t / 2).
    $ <eq_off_resonance>
    For the same-state transitions $ketbra(i, i)$ to $ketbra(i, i)$ the on- and off-resonance transitions are equal to
    $
        bra(i) cal(T)_"on" (ketbra(i, i)) ket(i) &= - sum_(j != i) bra(j) cal(T)_"on" (ketbra(i, i)) ket(j) \
        " and " bra(i) cal(T)_"off" (ketbra(i, i)) ket(i) &= - sum_(j != i) bra(j) cal(T)_"off" (ketbra(i, i)) ket(j).
    $ <eq_same_state_transition_resonances>
] <def_transition>
We will now use these definitions to prove @thm_tsp_second_order_expansion.
#proof([of @thm_tsp_second_order_expansion])[
    The bulk of this proof will be based on straightforward reductions from @eq_el_gigante. To start we will first show that no off-diagonal elements are introduced to the density matrix. By taking the $(j,k)$ matrix element of the output from @eq_el_gigante we see
    $
        bra(j) cal(T)(ketbra(i, i))ket(k) &= sum_(l, m) e^(-beta lambda_E (m)) / (1 + e^(-beta lambda_E (m))) bra(j\, l) alpha^2 / 2 EE_G [ partial^2 / (diff alpha^2) Phi_G (ketbra(i\, m, i\, m)) |_(alpha = 0) ] ket(k\, l) \
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
        &= tilde(alpha)^2 sum_(k, l) q(k) sinc^2 (frac(Delta(i, k | j, l) t, 2)) \
        &= tilde(alpha)^2 (q(0) sinc^2 (frac(Delta(i, 0 | j, 0) t, 2)) + q(0) sinc^2 (frac(Delta(i, 0 | j, 1) t, 2))) \
        & quad + tilde(alpha)^2 (q(1) sinc^2 (frac(Delta(i, 1 | j, 0) t, 2)) + q(1) sinc^2 (frac(Delta(i, 1 | j, 1) t, 2))) \
        &= tilde(alpha)^2 (q(0) sinc^2 (frac(Delta_S (i,j) t, 2)) + q(0) sinc^2 (frac((Delta_S (i,j) - gamma) t, 2))) \
        & quad + tilde(alpha)^2 (q(1) sinc^2 (frac((Delta_S (i, j) + gamma) t, 2)) + q(1) sinc^2 (frac(Delta_S (i,j) t, 2))),
    $
    where we see that combining the terms with $sinc^2 (Delta_S(i,j) t / 2)$, as $q(0) + q(1) = 1$, we immediately get @eq_transition_terms_total.

    To classify these terms as on-resonance or off-resonance we will focus on the argument to the sinc function, which is of the form $Delta_S(i,j) t / 2$ or $(Delta_S(i,j) plus.minus gamma) t / 2$. The idea is that we will take $t$ large enough so that only the energy differences that are less than $delta_"min"$, as defined in @eq_delta_min_def, will be non-negligible. Clearly the term $tilde(alpha)^2 sinc^2 ( frac(Delta_S (i,j)t, 2) )$ will always be off-resonance, as $delta_"min" <= Delta_S (i,j)$.

    Now we have three terms to classify as either on-resonance or off-resonance, we refer to each term by their argument to the $sinc$ function. The first we can categorically declare as being off-resonance is the $Delta_S(i,j)$ term. By [??] we know $sinc^2(Delta_S(i,j) t / 2) <= 4 / (delta_"min"^2 t^2)$, which we will make arbitrarily small in later sections. The other two can only be classified as on or off resonance depending if $Delta_S(i,j)$ is positive or negative. If $i > j$ then we know that $Delta_S(i,j) >= 0$ and therefore $sinc^2((Delta_S(i,j) - gamma)t / 2)$ term can be close to 1 if $gamma approx Delta_S(i,j)$, which also shows the $Delta_S(i,j) + gamma$ term is off-resonance for all $gamma$. We say that the $Delta_S(i,j) - gamma$ term in this scenario is on-resonance if $|Delta_S(i,j) - gamma| <= delta_"min"$. This classification is best described symbolically as
    $
        i > j "and" |Delta(i, j) - gamma| <= delta_("min") ==>
        bra(j) cal(T)_"on" (ketbra(i, i)) ket(j) =
        tilde(alpha)^2 q(0) "sinc"^2 ( (Delta_S (i,j) - gamma)t / 2 ).
    $ <tmp_tsp_1>
    The $q(0)$ prefactor indicates that the ancilla started in it's low energy state and since $sinc^2$ is symmetric we can write the argument as $gamma Delta_S (i,j)$ which can be remembered as the ancilla gaining $gamma$ amount of energy and the system losing $Delta_S (i,j)$. In this scenario the $Delta_S (i,j) + gamma$ term is therefore put in the off-resonance map
    $
        i > j "and" |Delta_S (i,j) - gamma| <= delta_min \
        ==>
        bra(j) cal(T)_"off" (ketbra(i, i)) ket(j) =
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
        bra(i) cal(T)(ket(i)bra(i)) ket(i) &= - tilde(alpha)^2 sum_(k != i) bra(k) cal(T)(ketbra(i, i)) ket(k) \
        &= - tilde(alpha)^2 sum_(k != i) bra(k) (cal(T)_"on" (ketbra(i, i)) + cal(T)_"off" (ket(i)bra(i))) ket(k) \
        &=: bra(i) cal(T)_"on" (ketbra(i, i)) ket(i) + bra(i) cal(T)_"off" (ket(i)bra(i)) ket(i)
    $

    By plugging in @tmp_tsp_1 and @tmp_tsp_2 to the above we are done with the self-transition terms.
]

=== Markovian Dynamics and Error Terms <sec_tsp_markovian_dynamics>

Now that we have fully computed the significant contributors to the output of our channel $Phi$, we move on to characterize the behavior of the channel as a Markov chain with noise.
A Markov chain is a random process that involves a walker transitioning to vertices on a graph wherein the probability of transition does not depend on the history of the walker. Specifically, in this context we view the vertices in this graph as the eigenstates of the Hamiltonian. The repeated interaction model because of the lack of coherences in the weak coupling limit can be interpreted as a Markov process over these eigenstates with transitions probabilities given by the above analysis.

Specifically, the Markov chain is dictated by the $Phi (rho; 0)$ and $cal(T)_"on"$ terms in the weak-coupling expansion, for $[rho, H_S] = 0$ we showed that $Phi (rho; 0) = id (\ho)$, so from now on we will specifically only deal with such density matrices and characterize the zeroth order term as an identity map. As for the Markov chain, we will use normal font to denote matrices, such as $I$ for the identity matrix and $T$ for the transition term added on. We use $e_i$ to denote the basis vector associated with the quantum state $ketbra(i, i)$ and $p$ to denote the probability vector for $rho$ associated with its eigenvalues.

#lemma([Quantum Dynamics to Classical Markov Chain])[
    Let $T$ be the matrix defined by
    $
        e_i^tpose T e_j := bra(i) cal(T)_"on" (ketbra(j, j)) ket(i).
    $
    The matrix $I + T$ is a column stochastic matrix and models the Markovian dynamics of our thermalizing channel up to $O(alpha^2 t^2)$,
    $
        bra(j) (id + cal(T)_"on")^(compose L) (ketbra(i, i)) ket(j) = e_j^tpose (I + T)^L e_i.
    $
    By linearity of $id + cal(T)_"on"$ this identity extends to any diagonal density matrix input $rho = sum_i p(i) ketbra(i, i)$.
] <lem_tsp_quantum_to_classical>
#proof()[We prove this inductively on $L$. The base case of $L = 1$ is trivial from the defintion of $T$
    $
        bra(j) (id + cal(T)_("on"))(ketbra(i, i)) ket(j) = delta_(i,j) + bra(j) cal(T)_("on")(ketbra(i, i)) ket(j) = e_j^tpose (I + T) e_i.
    $

    For the inductive step we will rely on the fact that there are no off-diagonal elements for diagonal inputs.
    $
        bra(j) cal(T)_("on") (ketbra(i, i)) ket(k) = delta_(j,k) bra(j) cal(T)_("on") (ketbra(i, i)) ket(j) ==> bra(j) cal(T)_("on")^(compose L) (ketbra(i, i)) ket(k) = delta_(j,k) bra(j) cal(T)_(on)^(compose L) (ketbra(i, i)) ket(j).
    $

    This is again by induction where the case $L = 1$ is proved in Theorem and the inductive step is
    $
        bra(j) cal(T)_(on)^(compose L) (ketbra(i, i)) ket(k) &= bra(j) cal(T)_(on) ( cal(T)_(on)^(compose L - 1) (ketbra(i, i)) ) ket(k) \
        &= sum_(m, n) bra(j) cal(T)_(on) ( ketbra(m, m) cal(T)_(on)^(compose L - 1)(ketbra(i, i)) ketbra(n, n) ) ket(k) \
        &= sum_(m, n) delta_(m,n) bra(m) cal(T)_(on)^(compose L - 1)(ketbra(i, i)) ket(m) bra(j) cal(T)_(on) ( ketbra(m, m) ) ket(k) \
        &= sum_(m) bra(m) cal(T)_(on)^(compose L - 1)(ketbra(i, i)) ket(m) delta_(j,k) bra(j) cal(T)_(on) ( ketbra(m, m) ) ket(j) \
        &= delta_(j,k) bra(j) cal(T)_(on)^(compose L) (ketbra(i, i)) ket(j).
    $

    This argument points the way towards how we will prove the inductive step in our stochastic conversion, starting with
    $
        bra(j) (identity + cal(T)_(on))^(compose L)(ketbra(i, i)) ket(j) &= bra(j) ( (identity + cal(T)_(on))^(compose L - 1) (ketbra(i, i)) + cal(T)_(on) compose (identity + cal(T)_(on))^(compose L - 1) (ketbra(i, i)) ) ket(j) \
        &= e_j^tpose (identity + T)^(L - 1) e_i + bra(j) cal(T)_(on) compose (identity + cal(T)_(on))^(compose L - 1) (ketbra(i, i)) ket(j).
    $<eq_tsp_matrix_reloaded1>

    We can use the inductive hypothesis on the term on the left and we now have to break down the $cal(T)_(on)$ term.
    $
        bra(j) cal(T)_(on) compose (identity + cal(T)_on)^(compose L - 1) (ketbra(i, i)) ket(j) &= sum_(m, n) bra(j) cal(T)_(on) ( ketbra(m, m) (identity + cal(T)_(on))^(compose L - 1) (ketbra(i, i)) ketbra(n, n) ) ket(j) \
        &= sum_(m) bra(j) cal(T)_(on) ( ketbra(m, m) ) ket(j) e_m^tpose (I + T)^(L - 1) e_i \
        &= sum_m e_j^tpose T e_m e_m^tpose (I + T)^(L -1) e_i \
        &= e_j^tpose T(I + T)^(L-1) e_i.
    $

    Substituting this into @eq_tsp_matrix_reloaded1 yields
    $ bra(j) (identity + cal(T)_(on))^(compose L)(ketbra(i, i)) ket(j) = e_j^tpose (I + T)^(L) e_i. $

    Our final step in the proof is to show that $I + T$ is column-stochastic. This is straightforward from our definition of $T$
    $ sum_i e_i^tpose (I + T) e_j = 1 + sum_i bra(i) cal(T)_(on)(ketbra(j, j)) ket(i). $

    Now we use the fact that $bra(j) cal(T)_(on)(ketbra(j, j)) ket(j) = - sum_(i != j) bra(i) cal(T)_(on)(ketbra(j, j)) ket(i)$ from @eq_same_state_transition_resonances to conclude that $I + T$ is column stochastic.
]

Since we will be effectively reducing our quantum dynamics to classical dynamics over the eigenbasis for $H_S$ we will need bounds on the convergence of Markov chains. This is a very deep area of research, with many decades of results, so we point interested readers to the comprehensive book by Levin and Peres @levin2017markov. As we will be dealing with non-reversible Markov chains we unfortunately cannot use the relatively well-developed theory for reversible Markov chains. Luckily, we will only need the following theorem due to Jerison.
#theorem([Jerison's Markov Relaxation Theorem @jerison2013general])[
    Let $M : bb(R)^(N) -> bb(R)^(N)$ be an ergodic Markov transition matrix acting on an $N$ dimensional state space with absolute spectral gap $lambda_star := 1 - max_(i > 1) |lambda_i (M)|$, where the eigenvalues of $M$ are ordered $1 = lambda_1 >= lambda_2 >= dots >= lambda_N >= -1$. Given this gap, if the number of steps $L$ in the Markov chain satisfies the following bound
    $
        L &>= (N) / (lambda_(star)) ( 2 log (1) / (lambda_(star)) + 4(1 + log 2) + (1) / (N) (2 log ( (1) / (epsilon) ) - 1) ) =: (N) / (lambda_star) J,
    $
    where $J$ is the collection of logarithmic and constant terms that we will typically ignore in asymptotic notation, then the resulting state $M^L vec(x)$ is $epsilon$ close to the fixed point
    $
        forall arrow(x) " s.t. " x_i >= 0 " and " sum_i x_i = 1, quad norm(arrow(pi) - M^L arrow(x))_1 <= epsilon.
    $
    We use $arrow(pi)$ to denote the unique eigenvector of eigenvalue 1 for $M$.
]<thm_markov_chain_bound>

#h(5mm) Now that we have an idea of how long it takes for our Markov chain to converge to the fixed points we need to show which states are actually fixed points. We demonstrate that for finite $beta$ any fixed point must satisfy a summation of detailed-balance terms. This fixed point is unique if the Markov chain is ergodic, which we do not argue in this lemma as an arbitrary thermalization channel $Phi$ may not be ergodic. For the ground state limit of $beta -> oo$ we show that the Markov matrix $I + T$ is upper triangular, which is crucial to our analysis of the spectral gap of the Markov chain in later results. We also demonstrate that the ground state is a fixed point in this limit nearly trivially.
#lemma([Markov Chain Fixed Points])[
    Let $T$ be the transition matrix with sum zero columns $sum_j e_j^tpose T e_i$ for all $i$, negative diagonal entries $e_i^tpose T e_i <= 0$, and off-diagonals smaller than 1 $e_j^tpose T e_i <= 1$ for $j != i$, associated with the on-resonance term $cal(T)_(on)$ of an arbitrary thermalizing channel $Phi$. A vector $arrow(p)$ is a fixed point of the Markovian dynamics $I + T$ if and only if it is in the kernel of $T$. This holds for finite $beta$ if the following is satisfied for all $j$
    $
        sum_(i != j) (e^(-beta lambda_S (i))) / (partfun_S (beta)) e_j^tpose T e_i - (e^(-beta lambda_S (j))) / (partfun_S (beta)) e_i^tpose T e_j = 0.
    $<eq_tsp_detailed_balance>
    In the $beta -> infinity$ limit the ground state $e_1$ is a fixed point and $T$ is upper triangular
    $
        lim_(beta -> infinity) (I + T) e_1 = e_1 " and " i > j ==> lim_(beta -> infinity) e_i^tpose T e_j = 0.
    $
]<lem_fixed_points>
#proof()[
    To show that the thermal state is the fixed point of the zero knowledge thermalizing channel we need to show that
    $T arrow(p)_(beta) = 0$ and that the Markov chain is ergodic. Ergodicity will be easy to prove so we focus on showing that $T arrow(p)_(beta) = 0$. This condition can be expressed as
    $
        e_j^tpose T arrow(p)_(beta) = sum_i (e^(-beta lambda_S (i))) / (partfun_S (beta)) e_j^tpose T e_i = 0.
    $<eq_tsp_thermal_state_tmp_1>
    We can make a quick substitution as we know the diagonal elements must equal the sum of the remainder of the column
    $
        e_i^tpose T e_i = - sum_(j != i) e_j^tpose T e_i,
    $
    which we can then pull out the $i = j$ term from the sum in @eq_tsp_thermal_state_tmp_1
    $
        e_j^tpose T arrow(p)_(beta) &= sum_(i != j) (e^(-beta lambda_S (i))) / (partfun_S (beta)) e_j^tpose T e_i - (e^(-beta lambda_S (j))) / (partfun_S (beta)) sum_(k != j) e_k^tpose T e_j,
    $
    which is 0 if and only if $arrow(p)_(beta)$ is a fixed point of $I + T$.

    We now show the $beta -> infinity$ case. We can show that $T$ is upper triangular using @thm_tsp_second_order_expansion which gives us the on-resonance transition amplitude. We assume $i < j$, which implies $Delta_S (i,j) <= 0$, and get
    $
        lim_(beta -> infinity) e_j^tpose T e_i &= lim_(beta -> infinity) bra(j) cal(T)_(on)(ketbra(i, i)) ket(j) \
        &= tilde(alpha)^2 lim_(beta -> infinity) [ (e^(-beta gamma)) / (1 + e^(-beta gamma)) II[ |Delta_S (i,j) + gamma| <= delta_(min)] sinc^2 ((Delta_S (i, j) + gamma) t / (2) ) ] \
        &= tilde(alpha)^2 II[ |Delta_S (i,j) + gamma| <= delta_(min)] sinc^2 ((Delta_S (i, j) + gamma) t / 2 ) lim_(beta -> infinity) (e^(-beta gamma)) / (1 + e^(-beta gamma)) \
        &= 0.
    $
    This further shows that the ground state is a fixed point, as every other eigenvector must have higher energy and therefore all on-resonance transitions _from_ the ground state must be 0
    $
        lim_(beta -> infinity) e_1^tpose T e_1 &= lim_(beta -> infinity) bra(1) cal(T)_(on)(ketbra(1, 1)) ket(1) = - sum_(j > 1) lim_(beta -> infinity) bra(j) cal(T)_(on)(ketbra(1, 1)) ket(j) = 0.
    $
    This then shows that the ground state is fixed
    $
        (I + T) e_1 = e_1,
    $
    and completes the proof.
]

#h(5mm) Using the decomposition from @thm_tsp_second_order_expansion and intermediate expressions in its proof we can now show why the off-resonance map $cal(T)_"off"$ is named "off-resonance"; even in the worst case scenario of choosing a bad value of $gamma$ such that all terms in $cal(T)$ end up in $cal(T)_"off"$ the trace norm of its output is always controllably small via $alpha$.

#corollary()[The induced trace norm of the off-resonance map $cal(T)_"off" (rho)$, for all density matrices $rho$ such that $[rho, H_S] = 0$ and $dim >= 2$, is upper bounded for all choices of the environment Hamiltonian $gamma$ by
    $
        norm(cal(T)_"off" (rho))_1 <= (8 alpha^2) / (delta_min^2).
    $
]<cor_tsp_t_off_norm>
#proof()[
    This result follows from applying bounds on the sinc function from @lem_sinc_poly_approx (given in @sec_appendix_haar) to the worst-case scenario off-resonance terms given in @eq_off_resonance.
    $
        i != j ==> abs(bra(j)cal(T)_("off")(ketbra(i, i))ket(j)) &<= tilde(alpha)^2 (4) / (delta_(min)^2 t^2) ( 1 + q(0) + q(1) ) = (8 alpha^2) / (delta_(min)^2(dim + 1)).
    $
    This allows us to bound the off-resonance self-transition term in @thm_tsp_second_order_expansion
    as
    $
        abs(bra(i)cal(T)_("off")(ketbra(i, i))ket(i)) &= abs(- sum_(j != i) bra(j) cal(T)_("off")(ketbra(i, i)) ket(j)) <= (dim_S - 1) (8 alpha^2) / (delta_(min)^2 (dim + 1)) <= (4 alpha^2) / (delta_(min)^2).
    $
    Now we can use this, along with our no off-diagonal output elements of $cal(T)$, to compute the trace norm of the off-resonance map
    $
        norm(cal(T)_("off")(rho))_1 &= sum_(j) abs(bra(j) cal(T)_("off")(rho) ket(j)) \
        &<= sum_(i, j) rho_(i,i) abs(bra(j)cal(T)_("off")(ketbra(i, i))ket(j)) \
        &= sum_(i) rho_(i,i) (sum_(j != i) |bra(j) cal(T)_("off")(ketbra(i, i)) ket(j)| + |bra(i) cal(T)_("off")(ketbra(i, i))ket(i)| ) \
        &<= sum_(i) rho_(i,i) ( (dim_S - 1) (8 alpha^2) / (delta_(min)^2(dim + 1)) + (4 alpha^2) / (delta_(min)^2) ) \
        &<= (8 alpha^2) / (delta_(min)^2).
    $
]

#h(5mm) The last result in this section that we will need is a bound on the trace norm of the remainder term, which we state in the following theorem.
#theorem([Remainder Bound])[
    Let $R_(Phi)(rho)$ be the remainder term for the second-order Taylor series expansion
    of the quantum channel $Phi$ acting on an input state $rho$ about $alpha=0$ defined in @eq_tsp_phi_taylor_series
    where the Schätten 1-norm of the remainder operator is bounded by
    $
        norm(R_(Phi) (rho ; alpha))_1 <= (16 sqrt(2)) / (sqrt(pi)) dim_S (alpha t)^3.
    $
]<thm_remainder_bound>
The proof of the remainder bound follows from the triangle inequality and remainder bounds on Taylor series and is given in @sec_appendix_haar.

I'm thinking of including a "template" theorem that can be used to simplify the 4 proofs contained in the following section. Let $rho_"fix"$ denote the unique fixed point for a channel $EE_gamma [id + cal(T)_"on"^((gamma))]$ and $tilde(lambda_star)$ the spectral gap of the scaled transition matrix, so $tilde(alpha)^2 tilde(lambda_star) = lambda_star$. Then we have
$
    norm(rho_"fix" - (EE_gamma Phi_gamma)^(compose L) (rho))_1 &<= norm(rho_"fix" - (id + EE_gamma cal(T)_"on"^((gamma)))^(compose L))_1 + L norm(EE_gamma cal(T)_"off"^((gamma)) + R_Phi)_1 \
    &<= norm(rho_"fix" - (id + EE_gamma cal(T)_"on"^((gamma)))^(compose L))_1 + L (norm(EE_gamma cal(T)_"off"^((gamma)))_1 + norm(R_Phi)_1) \
    &<= norm(rho_"fix" - (id + EE_gamma cal(T)_"on"^((gamma)))^(compose L))_1 + L ((8 alpha^2) / delta_min^2 + (16 sqrt(2)) / sqrt(pi) dim_S (alpha t)^3).
$
So now in order to balance these terms we can set $alpha = 1\/(dim_S delta_min^2 t^3 )$ and the expression on the right becomes $L alpha^2 / delta_min^2 (8 + 16 sqrt(2 / pi)).$ Now using Jerison's theorem we can argue that
$
    L >= dim^2 / (alpha^2 t^2 tilde(lambda_star)) J
$
is sufficient to guarantee that the distance to the fixed point is $tilde(O)(epsilon)$.
Now we note that the right hand side forces us to require $L alpha^2 / delta_min^2 in tilde(O)(epsilon)$ holds only if
$
    (dim^2) / (delta_min^2 t^2 tilde(lambda_star) ) in tilde(O)(epsilon)
$
can be satisfied if $t = dim / (delta_min sqrt(epsilon tilde(lambda_star)))$. and then we are done.

== Single Qubit and Truncated Harmonic Oscillator <sec_tsp_oscillator>

The first system we study is the qubit $hilb_S = CC^2$. This system is simple enough that we can explicitly write the dynamics as a $2 times 2$ transition matrix, which makes it easy to compute required simulation times and easy for the reader to follow. Although this system could be viewed as a warmup to the more general systems in Section \ref{sec:general_systems}, as the proof techniques are very similar, we remark that this system does have some unique properties. The biggest difference is that we do not assume any kind of belief distribution of the eigenvalue gap $Delta$ of the system. We only require that a window of width $2 sigma$ is known that contains $Delta$. We can then characterize the runtime in terms of $sigma$ and in addition to determining runtime we find it determines an upper bound on the $beta$ that can be prepared at low error.

The other unique phenomenon with the single qubit scenario is that the total simulation time needed is _independent_ of $beta$. Although this may seem incorrect, as most existing thermal state preparation algorithms tend to scale at least linearly with $beta$, this is in fact a property of the underlying Markov chain. The rate of convergence of the Markov chain is dictated by the spectral gap, which for this system is shown to be $tilde(alpha)^2$. The only aspect of the Markov chain that changes with $beta$ is what the fixed point is and the Markov Relaxation @thm_markov_chain_bound provides relaxation guarantees regardless of initial or final state.

#theorem()[
    Let $H_S$ be an arbitrary single qubit Hamiltonian with eigenvalue gap $Delta$, $rho$ any input state that commutes with $H_S$, and $L$ the number of interactions simulated. Given a window of width $2 sigma$ that is promised to contain $Delta$ and satisfies the inequality
    $
        sigma <= min {epsilon / (2 beta), Delta sqrt(epsilon / 2)},
    $
    then the following parameter choices
    $
        alpha &= 1 / (t^3(Delta + sigma)^2),\
        t &in 1 / sigma [sqrt(1- sqrt(1 - (2 sigma^2) / (epsilon Delta^2))), sqrt(1 + sqrt(1 - (2 sigma^2) / (epsilon Delta^2)))], \
        "and " L &= ceil(10 / (alpha^2 t^2(1 - sigma^2 t^2 \/2)) (2 log(5 / (alpha^2 t^2 sinc^2(|Delta - gamma| t \/2))) \ & #h(1cm) + 4 log(2 e) - 1 / 2 + log(2 / epsilon))),
    $
    are sufficient to guarantee thermalization of the form $norm(rho_S (beta) - Phi^(compose L) (rho))_1 in tilde(O)(epsilon)$. In the limit as $sigma -> 0$, the total simulation time required scales as
    $
        lim_(sigma -> 0) L dot t in tilde(O) (1 / (Delta epsilon^(2.5))).
    $
] <thm_single_qubit>
#proof()[
    The proof will be structured into three parts. First, we will need a bound on how close the fixed point of the Markov chain is to the thermal state, because the fixed point is exactly the thermal state only when $gamma = Delta$ and our window of width $sigma$ is sufficiently small given our error budget. Second, once we have these bounds we then need to determine the number of interactions $L$ that will be necessary to reach the fixed point within trace distance $epsilon$. Lastly, we use this value of $L$ to bound the accumulative error from the off-resonance mapping $cal(T)_"off"$ and remainder term $R_Phi$.

    We start by breaking down the trace distance into three components, one for the fixed-point distance from the thermal state, one for the Markov dynamics distance to the fixed-point, and lastly the remainder terms
    $
        &norm(rho_S (beta; Delta) - Phi^(compose L)(rho))_1 \
        &<= norm(rho_S (beta; Delta) - rho_S (beta; gamma))_1 + norm(rho_S (beta; gamma) - Phi^(compose L)(rho))_1 \
        &<= norm(rho_S (beta; Delta) - rho_S (beta; gamma))_1 + norm(rho_S (beta; gamma) - (identity + cal(T)_("on"))^(compose L)(rho))_1 + norm((identity + cal(T)_("on"))^(compose L)(rho) - Phi^(compose L)(rho))_1 \
        &<= norm(rho_S (beta; Delta) - rho_S (beta; gamma))_1 + norm(rho_S (beta; gamma) - (identity + cal(T)_("on"))^(compose L)(rho))_1 + L(norm(cal(T)_("off")(rho))_1 + norm(R_Phi)_1).
    $ <eq_single_qubit_three_errors>
    We proceed with the leftmost term first. The trace distance can be computed explicitly for a single qubit state as
    $
        &norm(rho_S (beta; gamma) - rho_S (beta; Delta))_1 \
        &= abs(bra(1) rho_S (beta; gamma) ket(1) - bra(1) rho_S (beta; Delta) ket(1)) + abs(bra(2) rho_S (beta; gamma) ket(2) - bra(2) rho_S (beta; Delta) ket(2)) \
        &= abs(bra(1) rho_S (beta; gamma) ket(1) - bra(1) rho_S (beta; Delta) ket(1)) + abs(1 - bra(1) rho_S (beta; gamma) ket(1) - 1 + bra(1) rho_S (beta; Delta) ket(1)) \
        &= 2 abs(bra(1) rho_S (beta; gamma) ket(1) - bra(1) rho_S (beta; Delta) ket(1)).
    $ <eq_single_qubit_int_1>
    Now we expand $bracket(1, rho_S (beta; gamma), 1)$ about $gamma = Delta$
    $
        bracket(1, rho_S (beta; gamma), 1) = frac(1, 1 + e^(-beta gamma)) &= frac(1, 1 + e^(-beta Delta)) + (gamma - Delta) beta frac(1, 1 + e^(-beta gamma_star)) frac(e^(-beta gamma_star), 1 + e^(-beta gamma_star)) \
        &= bracket(1, rho_S (beta; Delta), 1) + (gamma - Delta) beta frac(1, 1 + e^(-beta gamma_star)) frac(e^(-beta gamma_star), 1 + e^(-beta gamma_star)),
    $
    where $gamma_star$ denotes the special value of $gamma$ that is guaranteed to make the above equation hold by Taylor's Remainder Theorem.
    Since the rightmost factors can be upper bounded by $frac(1, 1 + e^(-beta gamma_star)) frac(e^(-beta gamma_star), 1 + e^(-beta gamma_star)) <= 1$, this can be rearranged and plugged into @eq_single_qubit_int_1 to give the upper bound
    $
        norm(rho_S (beta, gamma) - rho_S (beta, Delta))_1 <= 2 beta abs(Delta - gamma).
    $
    Since we require this distance to be less than $epsilon$, we can upper bound $abs(Delta - gamma) <= sigma$ and require
    $
        sigma <= frac(epsilon, 2 beta).
    $ <eq_single_qubit_ineq_1>

    Now we move on to the second stage of the proof: computing the number of interactions needed to reach the fixed point of the Markov chain. As the Markov transition matrix is only $2 times 2$ we will compute it explicitly. To do so, we need the matrix elements for $T$, which can be computed using @thm_tsp_second_order_expansion
    $
        arrow(e)_1^tpose T arrow(e)_1 = bracket(1, cal(T)_("on")(ket(1)bra(1)), 1) &= -tilde(alpha)^2 frac(e^(-beta gamma), 1 + e^(-beta gamma)) sinc^2(frac((-Delta + gamma)t, 2))
        \
        arrow(e)_2^tpose T arrow(e)_1 = bracket(2, cal(T)_("on")(ket(1)bra(1)), 2) &= tilde(alpha)^2 frac(e^(-beta gamma), 1 + e^(-beta gamma)) sinc^2(frac((-Delta + gamma)t, 2)) \
        arrow(e)_1^tpose T arrow(e)_2 = bracket(1, cal(T)_("on")(ket(2)bra(2)), 1) &= tilde(alpha)^2 frac(1, 1 + e^(-beta gamma)) sinc^2(frac((Delta - gamma)t, 2)) \
        arrow(e)_2^tpose T arrow(e)_2 = bracket(2, cal(T)_("on")(ket(2)bra(2)), 2) &= -tilde(alpha)^2 frac(1, 1 + e^(-beta gamma)) sinc^2(frac((Delta - gamma)t, 2)).
    $ <eq_single_qubit_markov_4>
    This gives us the total Markov chain matrix as
    $
        I + T = mat(1, 0; 0, 1) + tilde(alpha)^2 sinc^2((Delta - gamma) t / 2) q0 mat(-e^(-beta gamma), 1; e^(-beta gamma), -1),
    $ <eq_markov_matrix_single_qubit_gamma>
    where it can be seen that the fixed point is
    $
        (I + T) arrow(p)_(beta, gamma) = arrow(p)_(beta, gamma) = q0 arrow(e)_1 + q1 arrow(e)_2.
    $

    #h(5mm) To show convergence we will need the spectral gap of @eq_markov_matrix_single_qubit_gamma, which is given as $lambda_star = tilde(alpha)^2 sinc^2((Delta - gamma) t / 2)$. Plugging this into the Markov Relaxation @thm_markov_chain_bound we can compute a lowe bound on the number of interactions needed
    $
        L >= (2 J) / (tilde(alpha)^2 sinc^2( abs(Delta - gamma) t / 2)),
    $ <eq_one_qubit_l_bound_1>
    where $J$ captures subleading logarithmic factors.

    ur next goal is to simplify these bounds so that we can propagate them to our final error requirements. We first use @lem_sinc_poly_approx to produce a bound on $sinc$ whenever $gamma$ is within our window and $|Delta - gamma| <= sigma$
    $
        sinc^2( |Delta - gamma| t / 2) >= 1 - (sigma^2 t^2) / 2,
    $
    provided that $t sigma <= sqrt(2)$ to make the bound meaningful. Recalling that the dimension of the system is 4, we can then create a new lower bound for $L$ by plugging this expression for sinc in to @eq_one_qubit_l_bound_1
    $
        L >= (10 J) / (alpha^2 t^2 (1 - sigma^2 t^2 \/ 2))
    $ <eq_one_qubit_l_bound_2>
    which is larger than our bound in @eq_one_qubit_l_bound_1. If we choose $L$ to be twice this bound we will for sure meet the Markov chain error requirements.

    The third stage of the proof utilizes the above bound on $L$ to bound the off-resonance and remainder terms. The magnitiude of the total off-resonance contributions are $L norm(cal(T)_"off")_1 <= 8 \/ Delta^2$, given by @cor_tsp_t_off_norm, and the remainder term is $L norm(R_Phi)_1 <= 32 sqrt(2\/pi) alpha t^3$ from @thm_remainder_bound. Setting $alpha = 1\/ (t^3 (Delta + sigma)^2) <= 1\/(t^3 Delta^2)$ allows us to make the following inequalities
    $
        L(norm(cal(T)_"off")_1 + norm(R_Phi)_1) &<= (20 J) / (t^2 (1 - sigma^2 t^2 \/ 2)) (8 / Delta^2 + 32 sqrt(2 / pi) alpha t^3) \
        &<= (20 J) / (t^2 Delta^2 (1 - sigma^2 t^2 \/ 2)) (8 + 32 sqrt(2 / pi)).
    $ <eq_one_qubit_l_bound_3>
    The last step is then to show that the above is $tilde(O)(epsilon)$. As $J$ contains only logarithmic factors, it is sufficient to show that there exists a $t$ such that $t^2(1 - sigma^2 t^2 \/2) <= epsilon$. Rearranging this expression reveals a quadratic in $t^2$that must satisfy the following
    $
        Delta^2 t^2 (1 - (sigma^2 t^2) / 2) - 1 / epsilon >= 0.
    $ <eq_single_qubit_tmp_3>
    The roots of this quadratic are
    $
        t^2 = 1 / sigma^2 (1 plus.minus sqrt(1 - (2 sigma^2) / (Delta^2 epsilon))),
    $
    meaning that if $t$ lies between these two roots then our bound in @eq_one_qubit_l_bound_3 is $tilde(O)(epsilon)$.

    The first observation to make about these roots is that we require $sigma <= Delta sqrt(epsilon \/2)$ in order to keep the roots real and not become complex. As $sigma -> 0$ we note that the larger root $1 / sigma^2 (1 + sqrt(1 - (2 sigma^2) / (Delta^2 epsilon)))$ approaches infinity and the smaller root approaches $1\/Delta^2 epsilon$. This means that @eq_one_qubit_l_bound_3 has valid solutions provided $sigma$ is sufficiently small. This means that we have successfully bounded all 3 error terms present in the original decomposition @eq_single_qubit_three_errors by $tilde(O)(epsilon)$. We have done so by setting the following parameters
    - $alpha = 1 / ((Delta + sigma)^2 t^3)$,
    - $t in 1 / sigma^2 [1 - sqrt(1 - (2 sigma^2) / (Delta^2 epsilon)) , 1 + sqrt(1 - (2 sigma^2) / (Delta^2 epsilon))]$,
    - and $L >= tilde(Omega) (1 / (alpha^2 t^2 (1 - sigma^2 t^2 \/2)))$.
    Substituting in derived expressions for these parameters is sufficient to yield the theorem statement.
]

=== Harmonic Oscillator <sec_tsp_harmonic_oscillator>
Now that we have explored the thermalization channel completely for the single qubit case we turn our attention to a more complicated system: a truncated harmonic oscillator. For this scenario we will assume that the oscillator gap, $Delta$, is known. This is mostly to simplify proofs of ergodicity and should not be an issue in practice, as evidenced by later theorems that show thermalization without eigenvalue knowledge. The reason behind this proof requirement is that by tuning $gamma$ to be the spectral gap we can create a ``ladder" transition matrix in which states can move one level up or down. The proof of ergodicity relies on this ladder. Once we remove knowledge of $Delta$ if $gamma$ has some probability of being close to $2 Delta$ this special ladder structure is destroyed. To avoid this annoyance and focus on the special structure granted by the harmonic oscillator we will assume $gamma = Delta$.

This system also represents a transition from the single qubit to more general settings discussed later as the guarantees on total simulation time as a function of $beta$ are similar. For the harmonic oscillator we are only able to bound the spectral gap in the ground state limit as $beta -> infinity$, meaning that the convergence time for finite $beta$ has to be characterized in terms of the spectral gap of the Markov chain. For infinite $beta$ we are able to compute the spectral gap exactly, as the Markov transition matrix is upper triangular. The following theorem introduces this technique in a straightforward setting before it is used later for more complicated transition matrices.

#theorem([Harmonic Oscillator])[
    Let $H_S$ denote a truncated harmonic oscillator with $dim_S$ energy levels that are separated by $Delta$, giving $lambda_S (k) = k Delta$ for $1 <= k <= dim_S$, let $gamma$ be chosen to match the eigenvalue gap $gamma = Delta$ exactly, and let $rho$ be any input state that commutes with $H_S$. Setting the following parameters for the thermalizing channel $Phi$
    $
        alpha = (epsilon^(1.5) tilde(lambda_star)(beta)^(1.5) Delta) / (dim_S^4), t = dim_S (Delta sqrt(epsilon tilde(lambda_star)(beta)))," and " L in tilde(O)(dim_S^2 / (alpha^2 t^2 tilde(lambda_star)(beta))),
    $
    where $tilde(lambda_star)(beta)$ is the spectral gap of the scaled transition matrix $T \/ tilde(alpha)^2$, is sufficient for thermaliziation for arbitrary $beta$ as
    $
        norm(rho_S (beta) - Phi^(compose L) (rho))_1 in tilde(O)(epsilon).
    $
    This gives the total simulation time required as
    $
        L dot t in tilde(O)(dim_S^9 / (Delta epsilon^(2.5) tilde(lambda_star)(beta)^(2.5))).
    $
    In the limit $beta -> oo$ the above settings for $alpha, t$, and $L$ are valid for preparing the ground state with the spectral gap of the scaled transition matrix is further given by
    $
        lim_(beta -> oo) tilde(lambda_star)_beta = 1.
    $
] <thm_harmonic_oscillator>
#proof()[
    We first show that the thermal state is the unique fixed point for finite $beta$. This will be done by computing the nonzero on-resonance transitions and plugging in to @lem_fixed_points. As $gamma = Delta$, $Delta_S(i,j) = (i - j) Delta$, and $delta_min = Delta$, we can deduce that the on-resonance transitions will only be nonzero for adjacent states $ketbra(i, i)$ and $ketbra(i plus.minus 1, i plus.minus 1)$. This can be seen explicitly for $i != j$ by evaluating the transition elements given by @def_transition.
    $
        bra(j) cal(T)_"on" (ketbra(i, i)) ket(j)
        &= tilde(alpha)^2 q0 II[ abs(Delta_S (i,j) - gamma) <= delta_min] sinc^2((Delta_S (i,j) - gamma)t / 2) \
        &" " + tilde(alpha)^2 q1 II[abs(Delta_S (i,j) + gamma) <= delta_min] sinc^2((Delta_S (i,j) + gamma) t / 2)\
        &= tilde(alpha)^2 q(0) II[j = i - 1] sinc^2((Delta_S (i,j) - gamma)t / 2) \
        &" " + tilde(alpha)^2 q(1) II[j = i + 1] sinc^2((Delta_S (i,j) + gamma)t / 2) \
        &= tilde(alpha)^2 (q(0) II[j = i - 1] + q(1) II[j = i + 1])
    $ <eq_harmonic_oscillator_t_matrix>

    We now plug this expression into @eq_tsp_detailed_balance of @lem_fixed_points and use the fact that $Delta_S (i,i+1) = Delta$ for the harmonic oscillator
    $
        &sum_(i != j) frac(e^(-beta lambda_S (i)), cal(Z)_S (beta)) bracket(j, cal(T)_("on") (ket(i)bra(i)), j) - frac(e^(-beta lambda_S (j)), cal(Z)_S (beta)) bracket(i, TT_("on")(ket(j)bra(j)), i) \
        &= tilde(alpha)^2 e^(-beta lambda_S (j)) / (cal(Z)_S (beta) ) sum_(i != j) [II[j = i - 1] (q(0) e^(-beta Delta_S (i, j)) - q(1) ) + II[j = i + 1] (q(1) e^(-beta Delta_S (i, j)) - q(0) ) ] \
        &= tilde(alpha)^2 frac(e^(-beta lambda_S (j)), cal(Z)_S (beta)) ((e^(-beta Delta) q(0) - q(1)) + (e^(+beta Delta) q(1) - q(1))) \
        &= tilde(alpha)^2 frac(e^(-beta lambda_S (j)), cal(Z)_S (beta)) ((e^(-beta Delta) frac(1, 1 + e^(-beta gamma)) - frac(e^(-beta gamma), 1 + e^(-beta gamma))) + (e^(+beta Delta) frac(e^(-beta gamma), 1 + e^(-beta gamma)) - frac(1, 1 + e^(-beta gamma)))) \
        &= 0,
    $
    where the final equality comes from setting $gamma = Delta$. By @lem_fixed_points this is sufficient for $rho_S (beta)$ to be a fixed point of $id + cal(T)_"on"$.

    To show that $rho_S (beta)$ is the unique fixed point of the Markov chain it suffices to show that the walk is ergodic. This means that we need to show that the walk can generate transitions between any two sites, or in other words, the hitting time for any two states $i != j$ is nonzero. We prove this by induction on $i - j$ first for $i > j$. For $i = j + 1$ we have
    $
        bra(j) (id + cal(T)_"on")(ketbra(i, i)) ket(j) = bra(j) (id + cal(T)_"on")(ketbra(j + 1, j + 1)) ket(j) = tilde(alpha)^2 q(0),
    $
    which is nonzero and therefore the base case holds. Assuming $i = j + n$ holds we show that the hitting time for $i = j + n + 1$ is nonzero. Let $p$ denote the probability of transitioning from $j$ to $j + n$ after $n$ applications of $id + cal(T)_"on"$. We show that the probability of transitioning to $j + n + 1$ is nonzero in a few steps starting with the reduction
    $
        &bra(j) (id + cal(T)_"on")^(compose n + 1) (ketbra(j + n + 1, j + n + 1)) ket(j) \
        &= sum_(k_1, k_2) bra(j) (id + cal(T)_"on")^(compose n ) compose (ketbra(k_1, k_1) (id + cal(T)_"on")(ketbra(j + n + 1, j + n + 1)) ketbra(k_2, k_2)) ket(j).
    $
    We then can set $k_1 = k_2 = k$ as we know that $id + cal(T)_"on"$ does not add coherences, meaning it maps diagonal operators ($ketbra(j + n + 1, j + n + 1)$) to diagonal operators ($ketbra(k, k)$).
    We can then use the fact that one application of $id + cal(T)_"on"$ can map $j + n + 1$ to $j + n$ and take only that term out of the sum over $k$
    $
        &sum_k bra(j) (id + cal(T)_"on")^(compose n) compose (ketbra(k, k) (id + cal(T)_"on")(ketbra(j + n + 1, j + n + 1)) ketbra(k, k)) ket(j) \
        &>= bra(j) (id + cal(T)_"on")^(compose n) compose (braket(j + n, j + n) (id + cal(T)_"on")(ketbra(j + n + 1, j + n + 1)) ketbra(j + n, j + n)) ket(j) \
        &= tilde(alpha)^2 q(0) bra(j) (id + cal(T)_"on")^(compose n) (ketbra(j + n, j + n)) ket(j) \
        &= tilde(alpha)^2 q(0) p,
    $
    which is clearly greater than 0. To prove the case where $i < j$ the same inductive argument above can be repeated but this time factors of $q(1)$ accumulate as opposed to $q(0)$.

    ow that we have shown that the thermal state is the fixed point we would like to bound the total simulation time needed. To do so we first decompose our error into two parts, a Markov chain error and an off-resonance and remainder error
    $
        norm(rho_S (beta) - Phi^(compose L)(rho))_1 <= norm(rho_S (beta) - (id + cal(T)_"on" )^(compose L) (rho))_1 + L (norm(cal(T)_"off")_1 + norm(R_Phi)_1)
    $<eq_harmonic_oscillator_error_breakdown>
    We first bound the number of interactions, $L$, needed for the output of the Markov chain to be $epsilon$ close to the fixed point and then use this bound on $L$ to upper bound the off-resonance and remainder error. Unfortunately in the finite $beta$ scenario we are unable to determine the spectral gap of $T$, the entries of which are given in Eq. \eqref{eq:harmonic_oscillator_t_matrix}. The spectral gap of $T$ is necessary to use Jerison's Markov Relaxation @thm_markov_chain_bound which poses a problem for our understanding of the evolution time needed. Instead, we will pull out the overall factor of $tilde(alpha)^2$ and let $tilde(lambda_star) (beta)$ denote the spectral gap of $T \/ tilde(alpha)^2$. This then allows us to use @thm_markov_chain_bound but we will have to leave the number of interactions required in terms of $tilde(lambda_star)(beta)$.

    @thm_markov_chain_bound tells us that requiring
    $
        L >= dim_S / (tilde(alpha)^2 tilde(lambda_star) (beta)) J in tilde(O)(dim_S^2 / (alpha^2 t^2 tilde(lambda_star) (beta)))
    $

    is sufficient for the total variational distance between the stationary distribution to be $epsilon$-small, in other words $norm(rho_S (beta) - (id + cal(T)_"on")^(compose L)(rho))_1 in tilde(O)(epsilon)$. Now we use this expression for $L$ to bound the off-resonance and remainder errors. To do so we first want to asymptotically bound the two contributions, which can be found in @cor_tsp_t_off_norm and @thm_remainder_bound. The sum of the two errors is given by
    $
        norm(cal(T)_"off")_1 + norm(R_Phi)_1 <= (8 alpha^2) / (Delta^2) + 16 sqrt(pi / 2) dim_S (alpha t)^3.
    $
    by setting $alpha = 1\/(dim_S Delta^2 t^3)$ we can simplify the above as
    $
        norm(cal(T)_"off")_1 + norm(R_Phi)_1 <= (alpha^2) / (Delta^2) (8 + 16 sqrt(pi / 2)).
    $
    Using the sub-additivity property of the trace distance the total error scales as
    $
        L(norm(cal(T)_"off")_1 + norm(R_Phi)_1) &<= (dim_S^2 alpha^2 J) / (alpha^2 t^2 tilde(lambda_star)(beta) Delta^2) (8 + 16 sqrt(pi / 2)) \
        &in tilde(O)(dim_S^2 / (t^2 Delta^2 tilde(lambda_star) (beta))).
    $
    We can make this $tilde(O)(epsilon)$ by setting $t = dim_S / (Delta sqrt(epsilon tilde(lambda_star) (beta)))$. This then gives the total simulation time as
    $
        L dot t in tilde(O)(dim_S^2 / (alpha^2 t tilde(lambda_star) (beta) )) = tilde(O)(dim_S^9 / ( epsilon^(2.5) tilde(lambda_star) (beta)^(3.5))).
    $

    Now that we have analyzed the finite $beta$ regime, we turn to the $beta -> oo$ limit. Our proof above fo the fixed points only works for finite $beta$, but @lem_fixed_points tells us that in the $beta -> oo$ limit the ground state is a fixed point. We will show it is the unique fixed point by directly computing the spectrum of $T$, which will be rather easy to do. @lem_fixed_points further tells us that as $beta -> oo$ the matrix $T$ is upper triangular, which means we can compute the spectrum simply by just computing the diagonal elements. We do so via @eq_harmonic_oscillator_t_matrix and @eq_env_state_def, which says that for $1 < i < dim_S$ we have
    $
        arrow(e)_i^tpose T arrow(e)_i &= bra(i) cal(T)_"on" (ketbra(i, i)) ket(i) \
        &= - sum_(j != i) bra(j) cal(T)_"on" (ketbra(i, i)) ket(j) \
        &= - tilde(alpha)^2 sum_(j != i) (q(0) II[j = i - 1] + q(1) II[j = i + 1]) \
        &= -tilde(alpha)^2 (q(0) + q(1)) \
        &= - tilde(alpha)^2,
    $
    where the summation is only nonzero for $j = i plus.minus 1$. For $i = 1$ we note that because $arrow(e)_1$ is a fixed point we have $arrow(e)_j^tpose T arrow(e)_1 = 0$, so the diagonal entry is 0. The computation for $i = dim_S$ is similar to the above but yields from @eq_env_state_def
    $
        lim_(beta -> oo) bra(dim_S) cal(T)_"on" (ketbra(dim_S, dim_S)) ket(dim_S) = -tilde(alpha)^2 lim_(beta -> oo) q(0) = -tilde(alpha)^2.
    $

    #h(5mm) This shows us that the zero temperature limit of the transition matrix $T$ is
    $
        lim_(beta -> oo) T = tilde(alpha)^2 mat(
            0, 1, , ;
            , -1, 1, , ;
            , , -1, , ;
            , , , dots.down, ;
            , , , , 1;
            , , , , -1
        ).
    $
    We can compute the spectrum via the characteristic polynomial $det(lambda I - T)$. This is because $T$ is upper triangular and the determinant we need to compute is
    $
        det(lambda I - lim_(beta -> oo) T) = mat(
            delim: "|",
            lambda, -tilde(alpha)^2, , ;
            , lambda + tilde(alpha)^2, -tilde(alpha)^2, , ;
            , , lambda + tilde(alpha)^2, , ;
            , , , dots.down, ;
            , , , , -tilde(alpha)^2;
            , , , , lambda + tilde(alpha)^2
        ).
    $
    The roots of the above characteristic polynomial gives the spectrum of $T$ as 0 and $-tilde(alpha)^2$ with multiplicity $dim_S - 1$. This not only gives the spectral gap of $tilde(alpha)^2$ but further shows that the ground state is the unique fixed point because 0 only has multiplicity 1. This shows that $lim_(beta -> oo) tilde(lambda_star)(beta) = 1$.

    We now can use this to repeat the simulation time bound arguments from the finite $beta$ case. The decomposition in @eq_harmonic_oscillator_error_breakdown is still valid and we can use the Markov Relaxation @thm_markov_chain_bound to bound
    $
        L >= dim_S / (tilde(alpha)^2 lim_(beta -> oo) tilde(lambda_star) (beta)) J in tilde(Theta) (dim_S^2 / (alpha^2 t^2)).
    $
    The arguments for the off-resonance and remainder error bounds are the exact same and tell us that it suffices to set
    $
        alpha = 1 / (dim_S Delta^2 t^3) " and " t = dim_S / (Delta sqrt(epsilon)).
    $
    This gives the total simulation time needed as
    $
        L dot t in tilde(O)(dim_S^9 / (epsilon^(2.5) Delta)).
    $
]

=== Numerics <sec_specific_numerics>
Now that we have rigorous bounds on each of the parameters $alpha, t$ and $L$ needed to prepare thermal states of simple systems, we turn to numerics to test these bounds. The first question we explore is how the total simulation time $L dot t$ behaves as a function of $alpha$ and $t$. After, we examine the dependence of the total simulation time on the inverse temperature $beta$ and we observe a Mpemba-like effect where we find higher temperature states can cool faster than lower temperature ones @auerbach1995supercooling. Finally, we demonstrate how our proof techniques could be leading to worse $epsilon$ scaling than appears numerically necessary. Throughout these experiments we have the same numeric method of starting with the maximally mixed state $rho_S (0)$ and performing a search on the minimal number of interactions needed for the mean trace distance over all samples to be less than the target $epsilon$. The number of samples is increased until the variance in the trace distance is less than an order of magnitude below the mean.

In @fig_tot_time_vs_single_time we explore the total simulation time needed to prepare a thermal state with $beta = 2.0$ and $epsilon = 0.05$ for a single qubit system. We plot the total simulation time $L dot t$ needed as a function of $t$ for various settings of $alpha$. We find that increasing both parameters tends to decrease the overall cost until a saturation point is reached, which is at a value of $t$ slightly larger $1 / alpha$. For a fixed value of $alpha$ this initial decrease in $L dot t$ is inverse with $t$, in agreement with our finding of $L in tilde(O)(t^(-2))$ in Eq. \eqref{eq:single_qubit_l_bound_2} for $sigma = 0$. However, this process of decreasing the cost by increasing $t$ can only scale so far and appears to run into a minimum number of interactions $L$ required to thermalize. After this saturation point $L dot t$ scales linearly with $t$, indicating that the number of interactions $L$ has reached a minimum.

Another major take away from @fig_tot_time_vs_single_time is that it demonstrates that our thermalizing channel is exceptionally robust beyond the weak-coupling expansion in which we can theoretically analyze it. The values of $alpha t$ used in the far right of the plot completely break our weak-coupling expansion, as we have values of $tilde(alpha)$ that reach up to 500. One interesting phenomenon that we do not have an explanation for is the ``clumping" of various settings of $alpha$ in the large $t$ limit. As $alpha t$ dictates the amount of time that the random interaction term $G$ is simulated for, it could be that once a minimum amount of randomness is added via this interaction it is no longer beneficial in causing transitions among system eigenstates.

#figure(
    image("tsp_numerics/single_qubit_tot_time_vs_t.svg"),
    caption: [Total simulation time for a single qubit system to reach within trace distance of $0.05$ of the thermal state for $beta = 2$ as a function of per-interaction simulation time $t$. The slope of the large $t$ asymptote is $approx$ 1.01.],
)<fig_tot_time_vs_single_time>

The next task we have is to examine the $beta$ dependence. For the harmonic oscillator @thm_harmonic_oscillator is helpful for giving an idea of the total simulation time for the ground state but we cannot extend it to finite $beta$ due to the special structure of the transition matrix in the $beta -> oo$ limit. Perturbation theory could possibly be used to extend the computation of the spectral gap to the low temperature regime, but even then it would break down for large temperature (small $beta$). For generic $beta$ the structure of the harmonic oscillator transition matrix is tridiagonal but it is not quite Toeplitz, as the main diagonals deviate in the upper left and bottom right corners. We could try to pull these deviations into a separate matrix and treat them as perturbations to a fully Toeplitz matrix, which we can then compute the spectrum of. The issue with this approach is that these deviations are on the order of $tilde(alpha)^2 q(0)$ and $tilde(alpha)^2 q(1)$, which are comparable to the eigenvalues of the unperturbed matrix.

In @fig_sho_total_time_vs_beta we are able to probe the total simulation time and spectral gap of the harmonic oscillator as a function of $beta$. We reveal a rather surprising Mpemba-like phenomenon where it takes longer for an infinite temperature initial state (the maximally mixed state) to cool to intermediate temperatures than low temperature states. The Mpemba effect @mpemba is a classical phenomenon related to the time needed to freeze hot water compared to room temperature water with mentions going all the way back to Aristotle. This phenomenon has been extended to quantum thermodynamics and observed in both theory @nickMpemba @mpembaExplanation and in recent experimental research @zhang2025mpembaObservation. Our observations are not only a further analytic observation, but we are able to provide a proposed mechanism that explains the behavior. It is clear that the distance of our initial state to the target thermal state $norm(rho_S (beta) - rho_S (infinity))_1$ increases monotonically with $beta$ but what is not obvious is that the spectral gap of the underlying Markov chain is \emph{also} increasing. As larger spectral gaps lead to quicker convergences this acts in an opposite way on the total simulation time. The end result is that for small $beta$ the increase in initial distance is stronger than the increase in the spectral gap and $L dot t$ increases. After some amount of $beta$ these forces flip and the spectral gap effects become stronger than the initial state distance increasing, leading to a reduction in $L dot t$. This phenomenon appears to become more pronounced as the dimension of the harmonic oscillator increases, as can be seem in the $dim_S = 10$ data. Two things remain unclear: the first is what parameters affect the position and height of the peak in total simulation time and the second is if this behavior is present in Hamiltonians with more complicated eigenvalue difference structure than the harmonic oscillator.

#figure(
    grid(
        columns: 2,
        row-gutter: 3mm,
        image("tsp_numerics/sho_total_time_vs_beta_dim_4.svg"), image("tsp_numerics/spec_gap_dim_4.svg"),
        "(a) dim = 4 truncated harmonic oscillator", [(b) $dim = 4$ spectral gap vs $beta$],
        grid.cell(colspan: 2, image("tsp_numerics/sho_total_time_vs_beta_dim_10.svg", width: 50%)), grid.cell(
            colspan: 2,
            [(c) $dim = 10$ truncated harmonic oscillator],
        )
    ),

    caption: [Demonstration of $beta$ dependence of the thermalizing channel $Phi$ for the truncated harmonic oscillator. The environment gap $gamma$ was tuned to match the system gap $Delta$ exactly. The minimal number of interactions was found by binary search over values of $L$ that have an average error of less than $epsilon = 0.05$ with 100 samples.],
) <fig_sho_total_time_vs_beta>

The analytic proofs given in @thm_single_qubit and @thm_harmonic_oscillator are entirely based on our weak-coupling expansion derived in @sec_tsp_weak_coupling. The high level picture of this expansion is that we have a remainder error that scales like $tilde(O)(alpha t)^3)$ and an off-resonance error that scales as $O(alpha^2)$. To balance these two terms we then set $alpha = O(1 / t^3)$. However, as seen in @fig_tot_time_vs_single_time our thermalization routine appears to be quite robust beyond this weak-coupling expansion, which could lead to significant improvements in runtime. In our derivation for the $O(alpha)$ and $O(alpha^2)$ terms we relied on our eigenvalues being I.I.D Gaussian variables, with the first and second order expressions containing factors with the first and second moments respectively of the Gaussian distribution. This would suggest that the third order term in a weak coupling expansion might also be 0, similarly to the first order term. This would lead to a supposed remainder error of $O(alpha^4 t^4)$, which after balancing with the off-resonance error would give $alpha = O(1 \/ t^2)$. If the number of interactions then scales like $O(1 \/ (alpha^2 t^2))$, which is consistent with the spectral gap of $cal(T)_"on"$ scaling as $O(alpha^2 t^2)$, then to make the total error of order $O(epsilon)$ we would require $t in tilde(O)(1\/ epsilon^(0.5))$ as in @thm_single_qubit and @thm_harmonic_oscillator. This conjecture then leads to a total simulation time of order $O(1\/epsilon^(1.5))$.



#figure(
    image("tsp_numerics/epsilon_fitting_4.svg"),
    caption: [Scaling of $L dot t$ to prepare a harmonic oscillator thermal state with $beta = dim_S = 4$ with respect to $1 \/ epsilon$ in a log-log plot. For each line in the plot we scaled $alpha$ by a constant value to make $tilde(alpha)^2 approx 0.05$ for the largest value of $epsilon$. Each of these slopes was obtained via a least squares fitting of a power-law to $L dot t$ and $1 \/ epsilon$ and are consistently larger by 0.25-0.27 compared to stated predictions.],
) <fig_epsilon_scaling>

An even further conjecture would be to keep $alpha dot t$ as a small constant, in this case we are essentially saying that the randomized dynamics $e^(i alpha t G)$ are beneficial and should not be thought of as some remainder error to be minimized. If the $alpha t$ constant is small enough then the dynamics will still be approximated by the Markov chain $cal(T)_"on"$. Our spectral gap will still scale as $O((alpha t)^2)$ and $t$ as $O(1 \/ epsilon^(0.5))$. This would lead to our total simulation time scaling as $O(1 \/ epsilon^(0.5))$. In @fig_epsilon_scaling we numerically explore these various scalings of $alpha$ for the harmonic oscillator with $beta = dim_S = 4$. Our first remark is that the $alpha = O(1\/ t^3)$ scaling as dictated by @thm_harmonic_oscillator is numerically supported.
Specifically, the theorem suggests that we should observe $O(1 \/ epsilon^(2.5))$ scaling for $L dot t$.
Our experiment suggesting $L dot t in O(1 \/ epsilon^(2.764))$ which is approximately consistent and deviations from this scaling may arise from the inclusion of data in the fit from outside of the weak coupling limit which is the only regime where we anticipate this scaling.

== Generic Systems <sec_tsp_generic_sys>

We now extend our thermalization techniques to arbitrary Hamiltonians with no degenerate eigenvalues. The first major difficulty that we run into is how to choose our environment gap $gamma$. If one does not have any knowledge whatsoever about where the eigenvalues of $H_S$ may lie then we are reduced to uniform guessing. In Section \ref{sec:zero_knowledge} we show that even in this scenario the thermal state is an approximate fixed point for finite $beta$ and the exact fixed point for ground states and we provide a bound on the total simulation time required. However, we show that this generality does come at a cost. If one has complete knowledge of the eigenvalue differences we show in Section \ref{sec:perfect_knowledge} that the total simulation time markedly decreases. Further, with complete knowledge the thermal state is an exact fixed point for all $beta$. Finally, in Section \ref{sec:general_numerics} we study these impacts on small Hydrogen chain systems and observe the quantitative effects of noise added to $gamma$.

The assumption on non-degenerate eigenvalues is required for fairly technical conditions. In the $beta -> oo$ limit for our proof of the spectral gap we use the fact that the transition matrix $T$ is upper triangular in @lem_fixed_points. This is where the non-degeneracy is required because degenerate eigenvalues always have a non-zero transition amplitude that scales as $tilde(O)(alpha^2)$ without a factor of sinc. This means within the degenerate subspace in the transition matrix there is a uniform block. This makes computing the transition matrix spectrum a little more complicated than necessary, so we avoid this issue by requiring no degeneracies. This restriction could likely be lifted through an intelligent choice of eigenbasis for the degenerate subspace, or through better spectrum calculations of the resulting transition matrix, but we leave such explorations for future work.

=== Zero Knowledge <sec_tsp_zero_knowledge>

We now move on to show how our channel performs if one has no knowledge about the eigenvalue differences $Delta_S(i,j)$ apart from a bound on the maximum value of these differences. This is represented by choosing $gamma$ uniformly from the interval $[0, 4 norm(H_S)]$, which technically constitutes an upper bound on the largest $Delta_S (i,j)$, but estimates of $norm(H_S)$ are often readily attainable from the specification of the Hamiltonian using the triangle inequality. We also assume that an input state that commutes with the Hamiltonian can be provided, the maximally mixed state is sufficient as would a random eigenstate yielded by the quantum phase estimation algorithm.

#theorem([Zero Knowledge Thermal State Prep])[
    Let $H_S$ be a Hermitian matrix of dimension $"dim"_S$ with no degenerate eigenvalues, $rho$ any input state that commutes with $H_S$, and $gamma$ a random variable distributed uniformly in the interval $[0, 4 norm(H_S)]$ and let $rho_"fix"$ denote the unique fixed point of the transition dynamics $identity + EE_gamma cal(T)_"on"^((gamma))$ where $cal(T)_"on"^((gamma))$ is the on-resonance transition matrix used above with the dependence on $gamma$ made explicit. The following statements then hold.

    + For finite $beta$ the thermal state is an approximate fixed point of the thermalizing channel $EE_gamma Phi_gamma$ with a deviation of $ norm(rho_S (beta) - EE_gamma Phi_gamma (rho_S (beta)))_1 <= alpha^2 t e^(beta delta_"min") norm(H_S)^(-1) pi + 8 (alpha^2) / (delta_"min") + 16 sqrt((pi) / (2)) "dim"_S (alpha t)^3. $
    + The parameter settings $ alpha = (delta_"min"^4 epsilon^3 tilde(lambda)_star (beta)^3) / ("dim"_S^7 norm(H_S)^3), "  " t = ("dim"_S^2 norm(H_S)) / (epsilon tilde(lambda)_star (beta) delta_"min"^2), text("  and ") L in tilde(O)(("dim"_S^14 norm(H_S)^6) / (epsilon^5 delta_"min"^6 tilde(lambda)_star (beta)^6)) $ are sufficient for any $beta in [0,infinity]$ and error tolerance $epsilon in (0,2]$ to guarantee $norm(rho_"fix" - (EE_gamma Phi_gamma)^(compose L)(rho))_1 in tilde(O)(epsilon)$. The total simulation time needed is therefore $ L dot t in tilde(O) (("dim"_S^16 norm(H_S)^7) / (delta_"min"^8 epsilon^6 tilde(lambda)_star (beta)^7)). $
    + The fixed point is the ground state In the $beta -> infinity$ limit and the spectral gap, $tilde(lambda)_star (beta)$, of the rescaled transition matrix $EE_gamma T_gamma dot ((2 norm(H_S) ("dim" + 1)) / (alpha^2 t))$ is lower bounded by a constant, giving the two limits $ lim_(beta -> infinity) rho_"fix" = ket(1) bra(1) text(" and ") lim_(beta -> infinity) tilde(lambda)_star (beta) = 2 integral_0^(-delta_"min" t / 2) "sinc"^2(u) d u >= 2.43. $
] <thm_zero_knowledge>
#proof()[
    We start by understanding the fixed points of $id + EE_gamma cal(T)_(on)^((gamma))$, conditions for the thermal state being fixed are given in @lem_fixed_points. As the condition boils down to a detailed balance like condition, we need to compute the off-diagonal transition elements first. Starting with $i > j$ we have from @def_transition that
    $
        &EE_gamma bra(j) cal(T)_"on"^((gamma))(ketbra(i, i))ket(j) \
        &= tilde(alpha)^2 EE_(gamma) (1) / (1 + e^(-beta gamma)) II[ |Delta_S (i,j) - gamma| <= delta_"min"] "sinc"^2((Delta_S (i,j) - gamma)t / 2) \
        \
        &+ tilde(alpha)^2 EE_(gamma) (e^(-beta gamma)) / (1 + e^(-beta gamma)) II[ |Delta_S (i,j) + gamma| <= delta_"min"] "sinc"^2((Delta_S (i,j) + gamma)t / 2) \
        &= tilde(alpha)^2 EE_(gamma) (1) / (1 + e^(-beta gamma)) II[ |Delta_S (i,j) - gamma| <= delta_"min"] "sinc"^2((Delta_S (i,j) - gamma)t / 2) \
        &= tilde(alpha)^2 (1) / (4 norm(H_S)) integral_0^(4 norm(H_S)) (1) / (1 + e^(-beta gamma)) II [ |Delta_S (i,j) - gamma| <= delta_"min"] "sinc"^2((Delta_S (i,j) - gamma)t / 2) d gamma \
        &= (tilde(alpha)^2) / (4 norm(H_S)) integral_(Delta_S (i,j) - delta_"min")^(Delta_S (i,j) + delta_"min") (1) / (1 + e^(-beta gamma)) "sinc"^2((Delta_S (i,j) - gamma)t / 2) d gamma \
        &= (tilde(alpha)^2) / (2 t norm(H_S)) integral_(-delta_"min" t / 2)^(delta_"min" t / 2) (1) / (1 + e^(-beta (Delta_S (i, j) - 2 u / t))) sinc^2(u) d u
    $ <eq_zero_knowledge_transition_1>
    The exact same calculation holds for $i < j$, which after repeating the steps that led to @eq_zero_knowledge_transition_1 we arrive at a similar result with a slightly different integrand
    $
        EE_gamma bra(j) cal(T)_"on"^((gamma)) (ketbra(i, i)) ket(j) = tilde(alpha)^2 / (2 t norm(H_S)) integral_(- delta_min t / 2)^(delta_min t / 2) e^(- beta (Delta_S (j, i) - 2 u / t)) / (1 + e^(- beta (Delta_S (j, i) - 2 u / t))) sinc^2 (u) d u,
    $ <eq_zero_knowledge_transition_2>
    as we pick up a factor of $q(1)$ as opposed to $q(0)$. Note that we have also shown that $EE_gamma T_gamma$ is ergodic, as there is a nonzero probability for any state $ketbra(i, i)$ to transition to any other state $ketbra(j, j)$ in one iteration on average over $gamma$.

    For finite $beta$ the condition for $rho_S (beta)$ being a fixed point is given in @eq_tsp_detailed_balance, repeated here as
    $
        sum_(j != i) e^(-beta lambda_S (i)) / (cal(Z)_S (beta)) arrow(e)_j^tpose EE_gamma T_gamma arrow(e)_i - e^(-beta lambda_S (j)) / (cal(Z)_S (beta)) arrow(e)_i^tpose EE_gamma T_gamma arrow(e)_j = 0,
    $
    for all $j$. We can plug in our calculation for the transition coefficients for summands with $i > j$ first
    $
        &(e^(-beta lambda_S (i))) / (partfun_S (beta)) arrow(e)_j^tpose EE_gamma T_gamma arrow(e)_i - (e^(-beta lambda_S (j))) / (partfun_S (beta)) arrow(e)_i^T EE_gamma T_gamma arrow(e)_j \
        &= (e^(-beta lambda_S (j))) / (partfun_S (beta)) (e^(-beta Delta_S (i,j)) EE_gamma bra(j) cal(T)_"on"^((gamma))(ket(i)bra(i))ket(j) - bra(i) cal(T)_"on"^((gamma))(ket(j)bra(j))ket(i)) \
        &= (e^(-beta lambda_S (j))) / (partfun_S (beta)) (tilde(alpha)^2) / (2 t norm(H_S)) e^(-beta Delta_S (i,j)) integral_(-delta_"min" t / 2)^(delta_"min" t / 2) (1 - e^(beta 2 u / t)) / (1 + e^(-beta(Delta_S (i,j) - 2u / t))) "sinc"^2(u) d u \
        &= (e^(-beta lambda_S (i))) / (partfun_S (beta)) (tilde(alpha)^2) / (2 t norm(H_S)) integral_(-delta_"min" t / 2)^(delta_"min" t / 2) (1 - e^(beta 2 u / t)) / (1 + e^(-beta(Delta_S (i,j) - 2u / t))) "sinc"^2(u) d u.
    $ <eq_zero_knowledge_tmp_1>
    For $i < j$ we have the very similar
    $
        &(e^(-beta lambda_S (i))) / (partfun_S (beta)) arrow(e)_j^tpose EE_gamma T_gamma arrow(e)_i - (e^(-beta lambda_S (j))) / (partfun_S (beta)) arrow(e)_i^tpose EE_gamma T_gamma arrow(e)_j \
        &= (e^(-beta lambda_S (j))) / (partfun_S (beta)) (tilde(alpha)^2) / (2 t norm(H_S)) integral_(-delta_"min" t / 2)^(delta_"min" t / 2) (e^(beta 2 u / t) - 1) / (1 + e^(-beta (Delta_S (j, i) - 2u / t))) "sinc"^2(u) d u.
    $ <eq_zero_knowledge_tmp_2>
    Unfortunately these integrals are not 0, which can be verified numerically, and it is unclear how to make the summation over $i \neq j$ equal to 0.

    Our work around this is that instead of showing that the thermal state is exactly the fixed point we can use these results to show that it is an approximate fixed point. There are a few ways we could proceed. The first way could be to compute a Taylor series for the integrand and isolate the limits in which the remainder goes to 0. Unfortunately due to the $sinc^2(u) = sin(u)^2 / u^2$ term this means that the overall scaling will go like $1 / t$, making the total expression independent of $t$. Instead the route we will take will be to upper bound the norm $norm(arrow(p)_beta - EE_gamma (I + T_gamma) arrow(p)_beta)_1 = norm(EE_gamma T_gamma arrow(p)_beta)_1$, as this norm is only 0 if and only if $arrow(p)_beta$ is a fixed point. We reduce this to computations we have already performed as
    $
        norm(EE_gamma T_gamma arrow(p)_beta)_1 &= sum_j abs(arrow(e)_j^tpose EE_gamma T_gamma arrow(p)_beta) \
        &=sum_j abs(sum_i e^(-beta lambda_S (i)) / (cal(Z)_S (beta)) arrow(e)_j^tpose EE_gamma T_gamma arrow(e)_i) \
        &= sum_j abs(sum_(i != j) e^(-beta lambda_S (i)) / (cal(Z)_S (beta)) arrow(e)_j^tpose EE_gamma T_gamma arrow(e)_i - e^(-beta lambda_S (j)) / (cal(Z)_S (beta)) arrow(e)_i^tpose EE_gamma T_gamma arrow(e)_j) .
    $
    This is essentially the derivation for the fixed point conditions derived in @lem_fixed_points. We now plug in @eq_zero_knowledge_tmp_1 and @eq_zero_knowledge_tmp_2 into the above and upper bound the integral as
    $
        &sum_j abs(sum_(i != j) (e^(-beta lambda_S (i))) / (partfun_S (beta)) arrow(e)_j^tpose EE_gamma T_gamma arrow(e)_i - e^(-beta lambda_S (j)) / (cal(Z)_S (beta)) arrow(e)_i^tpose EE_gamma T_gamma arrow(e)_j) \
        &<= sum_j abs(sum_(i < j) (e^(-beta lambda_S (i))) / (partfun_S (beta)) arrow(e)_j^tpose EE_gamma T_gamma arrow(e)_i - e^(-beta lambda_S (j)) / (cal(Z)_S (beta)) arrow(e)_i^tpose EE_gamma T_gamma arrow(e)_j) + sum_j abs(sum_(i > j) (e^(-beta lambda_S (i))) / (partfun_S (beta)) arrow(e)_j^tpose EE_gamma T_gamma arrow(e)_i - e^(-beta lambda_S (j)) / (cal(Z)_S (beta)) arrow(e)_i^tpose EE_gamma T_gamma arrow(e)_j) \
        &= (tilde(alpha)^2) / (2 t norm(H_S)) sum_j abs(sum_(i < j) (e^(-beta lambda_S (j))) / (partfun_S (beta)) integral_(-delta_"min" t / 2)^(delta_"min" t / 2) (e^(beta 2 u / t) - 1) / (1 + e^(-beta(Delta_S (j, i) - 2u / t))) "sinc"^2(u) d u) \
        &+ (tilde(alpha)^2) / (2 t norm(H_S)) sum_j abs(sum_(i > j) (e^(-beta lambda_S (i))) / (partfun_S (beta)) integral_(-delta_"min" t / 2)^(delta_"min" t / 2) (1 - e^(beta 2 u / t)) / (1 + e^(-beta(Delta_S (i,j) - 2u / t))) "sinc"^2(u) d u) \
        &<= (tilde(alpha)^2) / (2 t norm(H_S)) sum_j sum_(i < j) (e^(-beta lambda_S (j))) / (partfun_S (beta)) integral_(-delta_"min" t / 2)^(delta_"min" t / 2) abs((e^(beta 2 u / t) - 1) / (1 + e^(-beta(Delta_S (j, i) - 2u / t)))) "sinc"^2(u) d u \
        &+ (tilde(alpha)^2) / (2 t norm(H_S)) sum_j sum_(i > j) (e^(-beta lambda_S (i))) / (partfun_S (beta)) integral_(-delta_"min" t / 2)^(delta_"min" t / 2) abs((1 - e^(beta 2 u / t)) / (1 + e^(-beta(Delta_S (i,j) - 2u / t)))) "sinc"^2(u) d u \
        &<= (tilde(alpha)^2) / (2 t norm(H_S)) e^(beta delta_"min") integral_(-delta_"min"t / 2)^(delta_"min"t / 2) "sinc"^2(u) d u (sum_j (e^(-beta lambda_S (j))) / (partfun_S (beta)) sum_(i < j) 1 + sum_j sum_(i > j) (e^(-beta lambda_S (i))) / (partfun_S (beta))) \
        &<= (tilde(alpha)^2 "dim"_S) / (t norm(H_S)) e^(beta delta_"min") pi \
        &<= alpha^2 t e^(beta delta_"min") norm(H_S)^(-1) pi.
    $
    For this we can have $alpha t$, which represents the total simulation time multiplied by the strength of the random interaction $G$, be constant and still take $alpha -> 0$ to achieve arbitrarily small error.

    Now we turn to bounding the total simulation time. We will let $rho_"fix"$ denote the fixed point of the dynamics. As before, we break the error into two pieces
    $
        norm(rho_"fix" - (EE_gamma Phi_gamma)^(compose L) ( rho))_1 <= norm(rho_"fix" - (id + EE_gamma cal(T)_"on"^((gamma)))^(compose L) (rho))_1 + L (norm(cal(T)_"off")_1 + norm(R_Phi)_1).
    $
    Let $tilde(lambda_star) (beta)$ denote the spectral gap for the rescaled transition matrix $EE_gamma T_gamma dot ( (2 norm(H_S) (dim + 1)) / (alpha^2 t))$, as this is the dimensionful prefactor in front of the transitions derived in @eq_zero_knowledge_tmp_2 and @eq_zero_knowledge_tmp_2. Jerison's Markov Relaxation @thm_markov_chain_bound tells us that taking $L$ to satisfy
    $
        L >= dim_S / lambda_star J in tilde(O) ((dim_S^2 norm(H_S)) / (alpha^2 t tilde(lambda_star) (beta)))
    $ <eq_zero_knowledge_tmp_3>
    is sufficient to guarantee $norm(rho_"fix" - (id + EE_gamma cal(T)_"on"^((gamma)))^(compose L) (rho))_1 in tilde(O)(epsilon)$. Now we balance the off-resonance and remainder errors
    $
        norm(cal(T)_"off")_1 + norm(R_Phi)_1 <= (8 alpha^2) / (delta_min^2) + 16 sqrt(pi / 2) dim_S (alpha t)^3 = alpha^2 / delta_min^2 (8 + 16 sqrt(pi / 2)),
    $
    and we see that setting $alpha = 1\/ (dim_S delta_min^2 t^3)$ makes the parenthesis a constant. To bound the total off-resonance and remainder error we take the product
    $
        L(norm(cal(T)_"off")_1 + norm(R_Phi)_1) in tilde(O) ((dim_S^2 norm(H_S)) / (alpha^2 t tilde(lambda_star) (beta)) alpha^2 / delta_min^2 ) = tilde(O) ((dim_S^2 norm(H_S)) / (t delta_min^2 tilde(lambda_star) (beta)) ).
    $
    Observe that setting
    $
        t = (dim_S^2 norm(H_S)) / (epsilon delta_min^2 tilde(lambda_star) (beta))
    $
    is sufficient to make the above product $L(norm(cal(T)_"off")_1 + norm(R_Phi)_1) in tilde(O) (epsilon)$.

    We now turn to the $beta -> oo$ limit. For this we note that @lem_fixed_points guarantees that the ground state is a fixed point and that $EE_gamma T_gamma$ is upper triangular. We will show the ground state is unique by computing the spectrum of $EE_gamma T_gamma$. For this we take the $beta -> oo$ limit of the transitions in @eq_zero_knowledge_transition_1 and @eq_zero_knowledge_transition_2, which will give us the diagonal elements and then the spectrum. Starting with $i > j$ given in @eq_zero_knowledge_transition_1 we get
    $
        lim_(beta -> oo) EE_gamma bra(j) cal(T)_"on"^((gamma)) (ketbra(i, i)) ket(j) = tilde(alpha)^2 / (2 t norm(H_S)) integral_(-delta_min t / 2)^(delta_min t / 2) sinc^2 (u) d u
    $
    and for $i < j$ from @eq_zero_knowledge_transition_2 we have
    $
        lim_(beta -> oo) EE_gamma bra(j) cal(T)_"on"^((gamma)) (ketbra(i, i)) ket(j) = 0.
    $
    We denote the sinc integration above as
    $
        I_"sinc" (t) := integral_(-delta_min t / 2)^(delta_min t / 2) sinc^2 (u) d u,
    $
    and we will show later that this is constant for $dim_S >= 3$. Now these transitions allow us to compute the diagonal elements
    $
        lim_(beta -> oo) EE_gamma bra(i) cal(T)_"on"^((gamma)) (ketbra(i, i)) ket(i) &= - sum_(j != i) lim_(beta -> oo) EE_gamma bra(j) cal(T)_"on"^((gamma)) (ketbra(i, i)) ket(j) \
        &= - sum_(j < i) lim_(beta -> oo) EE_gamma bra(j) cal(T)_"on"^((gamma)) (ketbra(i, i)) ket(j) \
        &= - tilde(alpha)^2 / (2 t norm(H_S)) (i - 1) I_"sinc" (t).
    $
    This gives a spectrum for $EE_gamma T_gamma$ as 0 and $- tilde(alpha)^2 / (2 t norm(H_S)) (i - 1) I_(sinc)(t)$ for $i > 1$. This shows the ground state is the unique fixed point as 0 has multiplicity 1 in the spectrum. Further the spectral gap of the rescaled transition matrix $tilde(lambda_star)(beta)$ is then given by
    $
        lim_(beta -> oo) tilde(lambda_star)(beta) = I_"sinc" (t).
    $
    We can repeat the analysis for finding suitable values for $alpha, t,$ and $L$ to guarantee thermalization and we find that
    $
        alpha = 1 / (dim_S delta_min^2 t^3)", " t = (4 dim_S^2 norm(H_S)) / (epsilon delta_min^2)", and " L = tilde(O)((dim_S^2 norm(H_S)) / (alpha^2 t))
    $
    are sufficient to guarantee $norm(ketbra(1, 1) - (EE_gamma Phi_gamma)^(compose L) (rho))_1 in tilde(O)(epsilon)$. Substituting this into @eq_zero_knowledge_tmp_3 yields
    $
        L dot t in tilde(O)((dim_S^16 norm(H_S)^7) / (delta_min^8 epsilon^6 tilde(lambda_star) (beta)^7))
    $
    as stated in the second claim of the theorem.

    Our final task is to justify the third claim of the theorem, which involves showing that $I_(sinc)(t)$ as constant is valid. Using the choice of $t$ directly above
    $
        I_"sinc" (t) = integral_(- delta_min t / 2)^(delta_min t / 2 ) sinc^2(u) d u = 2 integral_0^(dim_S^2 4 norm(H_S) \/ (epsilon delta_min)) sinc^2 (u) d u.
    $ <eq_zero_knowledge_sinc_integral>
    Now we note that this integral is monotonic with respect to the upper limit of integration with a final value of $lim_(t -> oo) I_(sinc)(t) = pi$. We note that we can capture a significant amount of this integral by just requiring the upper limit to be greater than the first zero of sinc located at $pi / 2$, which is true if $epsilon <= (4 dim_S^2 norm(H_S)) / (pi delta_min)$. This value can be computed as $2 integral_0^(pi / 2) sinc^2(u) d u >= 2.43$. This can be guaranteed by noting that $epsilon$ can be at most 2, so the upper limit in @eq_zero_knowledge_sinc_integral is satisfied if
    $
        epsilon <= 2 <= 3^2 / pi <= dim_S^2 / pi <= (dim_S^2 4 norm(H_S)) / (pi delta_min),
    $
    as $delta_min <= 4 norm(H_S)$. This shows that for our choice of $t$ then $abs(I_"sinc" (t) - pi) <= 0.71$, rendering it asymptotically constant as claimed.
]

#h(5mm) There are a few points that need to be addressed with the above theorem. The first is that our proof of the approximate fixed point utilizes rather poor bounds, resulting in diverging behavior as $beta -> oo$. For finite $beta$ our bounds on the change in the thermal state scales as $e^(beta delta_min)$, which diverges as $beta$ goes to $oo$, but in this exact same limit we are able to show that the ground state is the \emph{exact} fixed point of the Markov chain. This clear divergence in approximation error is a result of loose bounds and could be a potential avenue for improvement. The second point we would like to address is the rather high asymptotic scaling. This is the byproduct of a few things, the most important of which is the introduction of a $1 / t$ in the reduction of the $sinc$ integral. This causes a downstream effect of increasing the degree of each asymptotic parameter. To improve this one would need some kind of knowledge of the eigenvalues to prevent a uniform integration of each $sinc$ term. We study the limiting case of this by assuming sample access to the exact eigenvalue differences $Delta_S (i,j)$ in @sec_tsp_perfect_knowledge and obtain much improved scaling. The second source of inflation in our asymptotic scaling could be our weak-coupling approach to studying the channel. As explored numerically in @sec_specific_numerics we find that using different $alpha$ scalings with respect to $t$ can greatly effect the $epsilon$ scaling of the total simulation time $L dot t$. A higher order analysis of this channel could lead to anytically better guarantees on the thermalization time required, even in this zero knowledge scenario.


=== Perfect Knowledge <sec_tsp_perfect_knowledge>

Oftentimes when studying a system some knowledge of the eigenvalue gaps may be present. Our goal in this section is to study the extreme case of this scenario where one has knowledge of the exact eigenvalue differences. This is unlikely to happen with realistic quantum materials but instead serves as an ideal scenario for our channel to benchmark the effects of eigenvalue knowledge. Further, we note for some computational tasks, such as amplitude amplification, the eigenvalues may be explicitly computable and the real task is to find the dominant eigenvectors. A more realistic model for studying the impacts of eigenvalue knowledge on the total simulation time might be to place Gaussians at each of the $Delta_S (i,j)$ values with some width $sigma$. This is the model we use for numeric investigations in @sec_tsp_generic_numerics, but we were unable to compute the total simulation time required analytically. We find that our model of perfect knowledge allows us to show a reduced total simulation time budget, with the ratio of zero knowledge to perfect knowledge scaling as $tilde(O) ( (norm(H_S)^7) / (delta_min^7 epsilon^(3.5) tilde(lambda_star)(beta)^(3.5) ) )$,
which gives an explicit worst-case simulation time bound for ground state preparation.

#theorem([Perfect Knowledge Thermal State Prep])[
    Let $H_S$ be a Hermitian matrix of dimension $"dim"_S$ with no degenerate eigenvalues, $rho$ any input state that commutes with $H_S$, and let $gamma$ be a random variable with distribution $"Pr"[gamma = Delta_S (i,j)] = (eta_Delta_S (i,j)) / (binom("dim"_S, 2))$ where $eta_(i,j)$ is the number of times a particular eigenvalue difference appears. For any $beta in [0,infinity]$ the thermal state can be prepared with controllable error

    $ norm(rho_S (beta) - (EE_gamma Phi_gamma)^(compose L)(rho))_1 in tilde(O)(epsilon) $

    with the following parameter settings

    $
        alpha &= (delta_"min" epsilon^(1.5) tilde(lambda)_star (beta)^(1.5)) / ("dim"_S^7), t = ("dim"_S^2) / (delta_"min" epsilon^(0.5) tilde(lambda)_star (beta)^(0.5)), text(" and ") L in tilde(O)(("dim"_S^14) / (epsilon^2 tilde(lambda)_star (beta)^3)),
    $

    where $tilde(lambda)_star (beta)$ is the spectral gap of the rescaled transition matrix $EE_gamma T_gamma dot (binom("dim"_S, 2)) / (tilde(alpha)^2)$.
    This gives the total simulation time required as

    $ L dot t in tilde(O)(("dim"_S^16) / (delta_"min" epsilon^(2.5) tilde(lambda)_star (beta)^(3.5))). $

    All of the above conditions hold in the ground state limit as $beta -> infinity$ and further we can compute a lower bound on the spectral gap of the rescaled transition matrix as

    $ lim_(beta -> infinity) tilde(lambda)_star (beta) = min_(i > 1) sum_(j < i) eta_Delta(i, j) >= 1. $
] <thm_perfect_knowledge>
#proof()[
    This proof structure is structurally similar to the proof of @thm_perfect_knowledge.
    To show that the thermal state is the fixed point we will need to compute transition factors of the form $EE_gamma bra(j)cal(T)_("on")^((gamma))(ket(i)bra(i))ket(j)$ for use in @lem_fixed_points. Using the on-resonance definition in @eq_on_resonance we have for $i > j$
    $
        &EE_gamma bra(j) cal(T)_("on")^((gamma))(ket(i)bra(i))ket(j) \
        &= tilde(alpha)^2 EE_(gamma) (1) / (1 + e^(-beta gamma)) II[ |Delta_S (i,j) - gamma| <= delta_"min"] "sinc"^2((Delta_S (i,j) - gamma)t / 2) \
        &+ tilde(alpha)^2 EE_(gamma) (e^(-beta gamma)) / (1 + e^(-beta gamma)) II[ |Delta_S (i,j) + gamma| <= delta_"min"] "sinc"^2((Delta_S (i,j) + gamma)t / 2) \
        &= tilde(alpha)^2 sum_(Delta_S (k,l)) "Pr"[gamma = Delta_S (k,l)] (II[ |Delta_S (i,j) - Delta_S (k,l)| <= delta_"min"]) / (1 + e^(-beta Delta_S (k,l))) "sinc"^2((Delta_S (i,j) - Delta_S (k,l))t / 2) \
        &= tilde(alpha)^2 (eta_Delta(i, j)) / (binom("dim"_S, 2)) (1) / (1 + e^(-beta Delta_S (i,j))).
    $
    $i < j$ can be computed similarly as
    $
        EE_gamma bra(j) cal(T)_("on")^((gamma))(ket(i)bra(i))ket(j) = tilde(alpha)^2 (eta_Delta(i, j)) / (binom("dim"_S, 2)) (e^(-beta Delta_S (k,l))) / (1 + e^(-beta Delta_S (k,l))).
    $

    This allows us to compute the detailed-balance like condition in @eq_tsp_detailed_balance for $i > j$
    $
        &( e^(-beta lambda_S (i))) / (partfun_S (beta)) EE_gamma bra(j) cal(T)_("on")^((gamma))(ket(i)bra(i)) ket(j) - (e^(-beta lambda_S (j))) / (partfun_S (beta)) bra(i) cal(T)_("on")^((gamma))(ket(j)bra(j)) ket(i) \
        &= (e^(-beta lambda_S (i))) / (partfun_S (beta)) tilde(alpha)^2 (eta_Delta(i, j)) / (binom("dim"_S, 2)) (1) / (1 + e^(-beta Delta_S (i,j))) - (e^(-beta lambda_S (j))) / (partfun_S (beta)) tilde(alpha)^2 (eta_Delta(i, j)) / (binom("dim"_S, 2)) (e^(-beta Delta_S (i,j))) / (1 + e^(-beta Delta_S (i,j))) \
        &= (tilde(alpha)^2) / (partfun_S (beta)) (eta_Delta(i, j)) / (binom("dim"_S, 2)) ((e^(-beta lambda_S (i))) / (1 + e^(-beta Delta_S (i,j))) - e^(-beta lambda_S (j)) (e^(-beta Delta_S (i,j))) / (1 + e^(-beta Delta_S (i,j)))) \
        &= 0.
    $

    For $i < j$ we can repeat the same steps to argue that detailed balance also holds in this case.

    $
        &(e^(-beta lambda_S (i))) / (partfun_S (beta)) EE_gamma bra(j) cal(T)_("on")^((gamma))(ket(i)bra(i)) ket(j) - (e^(-beta lambda_S (j))) / (partfun_S (beta)) bra(i) cal(T)_("on")^((gamma))(ket(j)bra(j)) ket(i) \
        &= (e^(-beta lambda_S (i))) / (partfun_S (beta)) tilde(alpha)^2 (eta_Delta(i, j)) / (binom("dim"_S, 2)) (e^(-beta Delta_S (j, i))) / (1 + e^(-beta Delta_S (j, i))) - (e^(-beta lambda_S (j))) / (partfun_S (beta)) tilde(alpha)^2 (eta_Delta(i, j)) / (binom("dim"_S, 2)) (1) / (1 + e^(-beta Delta_S (j, i))) \
        &= (tilde(alpha)^2) / (partfun_S (beta)) (eta_Delta(i, j)) / (binom("dim"_S, 2)) ((e^(-beta lambda_S (j))) / (1 + e^(-beta Delta_S (j, i))) - (e^(-beta lambda_S (j))) / (1 + e^(-beta Delta_S (j, i)))) \
        &= 0.
    $

    This is sufficient to show that the thermal state $rho_S (beta)$ is a fixed point via @lem_fixed_points. As we have also shown that the probability of transitioning from any state $ket(i)bra(i)$ to any other state $ket(j)bra(j)$ is nonzero this gives a nonzero expected hitting time for any pair of states. This implies the Markov chain is ergodic and that $rho_S (beta)$ is the _unique_ fixed point.

    Next we bound the total simulation time required. For reasons similar to the harmonic oscillator in @sec_tsp_harmonic_oscillator we are unable to compute the spectral gap of the Markov matrix. We start the analysis in a similar manner by using the decomposition
    $
        norm(rho_S (beta) - (EE_gamma Phi_gamma)^(compose L)(rho))_1 <= norm(rho_S (beta) - (EE_gamma identity + cal(T)_("on")^((gamma)))^(compose L)(rho))_1 + L(norm(TT_"off")_1 + norm(R_(Phi))_1).
    $

    We bound the Markov error via @thm_markov_chain_bound. This theorem guarantees that choosing $L$ to satisfy
    $
        L >= ("dim"_S binom("dim"_S, 2)) / (tilde(alpha)^2 tilde(lambda)_star (beta)) J in tilde(O)(("dim"_S^4) / (alpha^2 t^2 tilde(lambda)_star (beta))),
    $

    where $tilde(lambda)_star (beta)$ is the spectral gap of the rescaled transition matrix $EE_gamma T_gamma dot (binom("dim"_S, 2)) / (tilde(alpha)^2)$, is sufficient for $norm(rho_S (beta) - (EE_gamma identity + cal(T)_("on")^((gamma)))^(compose L)(rho))_1 in tilde(O)(epsilon)$. We now use this to bound the total off-resonance and remainder error after balancing the two contributions asymptotically

    $
        norm(cal(T)_"off")_1 + norm(R_(Phi))_1 <= (8alpha^2) / (delta_"min"^2) + 16 sqrt((pi) / (2)) "dim"_S (alpha t)^3 = (alpha^2) / (delta_"min"^2) (8 + 16 sqrt((pi) / (2)) alpha "dim"_S delta_"min"^2 t^3).
    $

    Setting $alpha = (1) / ("dim"_S delta_"min"^2 t^3)$ is sufficient to make the parenthesis a constant. Lastly to get the total error in $tilde(O)(epsilon)$ we multiply the above by the $L$ chosen before

    $
        L (norm(TT_"off")_1 + norm(R_(Phi))_1) in tilde(O)(("dim"_S^4) / (alpha^2 t^2 tilde(lambda)_star (beta)) (alpha^2) / (delta_"min"^2)) = tilde(O)(("dim"_S^4) / (t^2 delta_"min"^2 tilde(lambda)_star (beta))).
    $

    Choosing

    $ t = ("dim"_S^2) / (delta_"min" sqrt(epsilon tilde(lambda)_star (beta))) $

    is sufficient to guarantee $L (norm(TT_"off")_1 + norm(R_(Phi))_1) in tilde(O)(epsilon)$ and that the total error $norm(rho_S (beta) - (EE_gamma Phi_gamma)^(compose L)(rho))_1 in tilde(O)(epsilon)$. Combining the above results for $alpha, L$ and $t$ yields the theorem statement for finite $beta$.

    We now show how to calculate $tilde(lambda)_star (beta)$ in the $beta -> infinity$ limit. From @lem_fixed_points we know that $EE_gamma T_gamma$ will be upper triangular, implying again that we can compute the spectrum if we can compute the diagonal elements of the matrix. Using our computation of the off-diagonal elements from the proof of @thm_zero_knowledge we have for $i > 1$

    $
        lim_(beta -> infinity) EE_gamma bra(i) cal(T)_("on")^((gamma))(ket(i)bra(i)) ket(i) &= - lim_(beta -> infinity) sum_(j != i) bra(j) cal(T)_("on")^((gamma))(ket(i)bra(i)) ket(j) \
        &= - lim_(beta -> infinity) sum_(j < i) tilde(alpha)^2 (eta_Delta(i, j)) / (binom("dim"_S, 2)) (1) / (1 + e^(-beta Delta_S (i,j))) - lim_(beta -> infinity) sum_(j > i) tilde(alpha)^2 (eta_Delta(i, j)) / (binom("dim"_S, 2)) (e^(-beta Delta_S (j, i))) / (1 + e^(-beta Delta_S (j, i))) \
        &= - (tilde(alpha)^2) / (binom("dim"_S, 2)) sum_(j < i) eta_Delta(i, j).
    $

    For $i = 1$ as we know the ground state is fixed we have $lim_(beta -> infinity) EE_gamma bra(1) cal(T)_("on")^((gamma))(ket(1)bra(1)) ket(1) = 0$. This gives the spectrum of $EE_gamma T_gamma$ as 0 and $-(tilde(alpha)^2) / (binom("dim"_S, 2)) sum_(j < i) eta_Delta(i, j)$ for all $i > 1$. From this spectrum we can conclude that the ground state is the _unique_ fixed point as 0 has multiplicity 1 in the spectrum and further that the spectral gap can be bounded from below as

    $
        lim_(beta -> infinity) tilde(lambda)_star (beta) = min_(i > 1) sum_(j < i) eta_Delta(i, j) >= 1.
    $
]

The above theorem shows that if we sample our transitions strategically rather than randomly then we can achieve much faster convergence to the groundstates in our upper bounds. Importantly, the scaling of the total simulation time is also independent of the norm of $H_S$ in this case, whereas the time required by the zero knowledge case does. Unfortunately, the dimensional scaling of $dim_S^16$ is prohibitive for all but the smallest dimensional systems. This scaling is again likely loose because of a number of assumptions that we make above and also a result of our insistence that the channel always operate inside the regime of weak coupling. In contrast, we will see below that equilibration can be much faster if strong coupling is assumed. Finally, it is worth noting that although perfect knowledge is assumed, a cooling schedule is not used. By changing the distribution depending on the temperature of the Gibbs state it is possible that even better scaling may be achievable.

=== Numerics <sec_tsp_generic_numerics>

The analytic results developed in the previous two sections provide strong guarantees on the correctness of our routine for most quantum systems, however the bounds on the total simulation are fairly high degree polynomials in the parameters of interest. One crucial interpretation of the two different results is that knowledge of the eigenvalue differences of $H_S$ can lead to significantly better simulation time bounds, but this knowledge is not \emph{crucial} for thermalization. Another important takeaway is that we cannot bound the simulation time or number of interactions required for finite $beta$ as we cannot bound the spectral gap of the expected transition matrix $EE_gamma T_gamma$. The purpose of this section is to investigate these two theoretic takeaways numerically with small Hydrogen chain systems. These systems are some of the smallest chemical systems that still display some real-world chemical behavior, and as a result are typically used in many numeric benchmarks for quantum routines.

Our first experiment conducted is to study the effects of changing $alpha$ and $t$ on the trace distance error as a function of $L$. The theory developed in prior sections is very prescriptive; to reach a specific trace distance of $epsilon$ all of our theorems give a value of $alpha$, $L$, and $t$ that guarantee a distance of at most $tilde(O)(epsilon)$ but say nothing about what this convergence looks like. In Figure \ref{fig:h_chain_error} we study the effects of different choices of $alpha$ and $t$ on this convergence rate. To generate the Hamiltonians used in these experiments we created a small chain of equally spaced hydrogen nuclei with an STO-3G active space for the electrons. Hamiltonian creation was done with OpenFermion @mcclean2020openfermion and PySCF @pyscf. Once the Hamiltonians were generated, the distance to the thermal state $rho_S (beta)$ for each was tracked over $L = 5000$ interactions. For both Hydrogen 2 and Hydrogen 3 we chose $beta = 4$ for consistency, this gave a ground state overlap of around 0.56 for Hydrogen 2 and 0.26 for Hydrogen 3.
#todo[Need to replace this with actual plots.]
#figure(
    grid(
        columns: 2,
        row-gutter: 3mm,
        image("tsp_numerics/error_vs_interaction_h2_chain_1.svg"),
        image("tsp_numerics/error_vs_interaction_h3_chain_3.svg"),

        [(a)], [(b)],
    ),
    caption: [
        These plots show the distance to the target thermal state for Hydrogen 2 and Hydrogen 3 chains as the number of interactions $L$ increases. For both Hydrogen 2 and 3 we set $beta = 4.0$, which gives a ground state overlap of greater than 0.5 for Hydrogen 2 and 0.25 for Hydrogen 3. $gamma$ for both (a) and (b)) was generated by placing a Gaussian at the average energy $tr(H_S) / dim_S$ with a width of $norm(H_S) / 2$. We note that a variety of $tilde(alpha)^2$ values were chosen to demonstrate the faster convergence, but higher error, of strong coupling.
    ],
) <fig_h_chain_error>


There are a few key takeaways from Figure \ref{fig:h_chain_error}. The first is that we observe increasing $alpha$ and $t$ tend to increase the convergence rate, with $alpha$ visuallly appearing more important. At higher values of $alpha$ changes in $t$ appear to make less of an impact on the error. We also observe that our channel is seemingly robust beyond our weak-interaction analysis. For values of $tilde(alpha)^2$, the weak-interaction expansion parameter, we observe that values as high as $tilde(alpha)^2 approx 3$ can have rapid convergence to fairly low error floors. We note that these coupling values that go beyond weak-interaction seem to lead to faster convergence of the dynamics at the cost of a larger error floor. It remains an open question if dynamic choices of $alpha$ and $t$ could lead to better performance of the overall routine, one could use very large $tilde(alpha)$ initially to quickly thermalize with large error and then decrease $tilde(alpha)$ to fine-tune the final state.

The second observation we make is on the choice of the environment gaps $gamma$. For both Hydrogen 2 and 3 we selected $gamma$ randomly from a Gaussian with mean $tr(H_S) / dim_S$ and standard deviation of $norm(H_S) / 2$. This choice of $gamma$ is completely heuristic and was intended to have a large overlap with what the typical eigenvalue differences may look like with a large enough deviation to pick up potentially large differences. Although this heuristic works well enough to show convergence, it leads us to question if the error convergence or floors can be improved with better choices of $gamma$.

In @fig_h_chain_noise we demonstrate that better choices of $gamma$ do in fact reduce the total simulation time needed for thermalization. To generate the data we compute the number of interactions needed at a fixed coupling constant $alpha$ and a fixed time $t$ as a function of the noise added to our samples for $gamma$. We generate one sample of $gamma$ by first computing the eigenvalue spectrum of H2 or H3 exactly, then by choosing two non-equal eigenvalues, and finally sampling a Gaussian centered at the absolute value of the difference. The width of this Gaussian then serves as a proxy for the amount of knowledge one may have about the system's eigenvalues. We plot the total simulation time with respect to this width as it varies from 0 to the spectral width $max_i lambda_S (i) - min_j lambda_S (j)$. The results align well with our theoretic analysis: having knowledge of the eigenvalues of the system can be used to speed up the thermalization routine but if one does not have any knowledge at all the thermal state can still be prepared. It is an open question if the dependence of the total simulation time $L dot t$ on the noise level $sigma$ can be determined analytically.

#figure(
    grid(
        columns: 2,
        row-gutter: 3mm,
        image("tsp_numerics/h2_chain_with_noise_3.svg"), image("tsp_numerics/h3_chain_with_noise_3.svg"),
        [(a) Hydrogen 2], [(b) Hydrogen 3],
    ),
    caption: [
        In these plots the amount of total simulation time needed to prepare a $beta = 2.0$ thermal state with $alpha = 0.01$ and $t = 500$ is tracked as a function of the noise added to samples of $gamma$. A sample for $gamma$ is generated by choosing two non-equal eigenvalues from the system spectrum and adding a Gaussian random variable with standard deviation $sigma$. For each value of $sigma$ the resulting state needs to have an average trace distance of less than $0.05$ for 100 samples.
    ],
) <fig_h_chain_noise>

== Discussion <sec_tsp_discussion>

Thermal state preparation is likely to be a crucial preparation step for the simulation of quantum systems on digital quantum computers. We have presented a thermalization routine for this task that has an optimally minimal number of overhead ancilla qubits and compiles to remarkably simple circuits of time independent Hamiltonian evolution of the unprocessed system Hamiltonian, with no filtering or rejection steps and no Fourier weighted jump operators of Linbladians. Our routine is based on relatively recent classical Monte Carlo techniques, specifically Hamiltonian Monte Carlo @hoffman2011nouturnsampleradaptivelysetting and the end result bears striking resemblance to the Repeated Interactions framework in open quantum systems @prositto2025equilibrium. In the Hamiltonian Monte Carlo algorithm thermal states over a position coordinate $q$ is prepared by sampling momentum $p$ from the Boltzmann distribution for Gaussians $e^(-beta p^2 / 2m)$ followed by time evolution. Classical Hamiltonian dynamics is enough to couple the position and momentum, leading to the Boltzmann distribution over $q$ with enough time and samples. In the repeated interactions framework a quantum system interacts with many small environments, typically a single photon, that is repeated until the system thermalizes.

Our work extends these procedures to quantum algorithms. For Hamiltonian Monte Carlo, instead of adding in momentum variables we add in a single ancilla qubit to serve as our extra state space. We do not have the luxury of classical Hamiltonian dynamics that couples these two spaces or registers, so we add in a randomized interaction term to the Hamiltonian. After simulating the time dynamics of this system-ancilla pair and repeating multiple times we are able to thermalize the system to the same $beta$. On the other hand, the Repeated Interactions framework typically is concerned with thermodynamic limits, such as infinite time or interactions, and specific system-interactions pairs. As our procedure is intended to be used as a subroutine for quantum computers our techniques work for arbitrary, non-degenerate, Hamiltonians and purposefully use randomized interactions as opposed to a fixed interaction model.

One benefit of our thermalization procedure is that it can be compiled all the way down to the circuit level with minimal overhead in complexity. The elements of the channel that need to be compiled are: the ancilla qubit state preparation, the initial state preparation for the system, the time evolution of $H + alpha G$, and the partial trace. The ancilla state preparation can be done starting from the ground state with a Pauli-$X$ rotation. The initial system state can be any state that commutes with $H_S$, the two leading choices are the maximally mixed state or a Haar random pure state. The maximally mixed state can be prepared using just CNOT gates at a cost of higher qubit count and a Haar random pure state can be prepared with no additional qubit overhead by applying a deep enough random circuit @choi2023preparing.

The time evolution of $H + alpha G$ can be broken down into simulation of $H_S$ and then $alpha G$ by one application of Trotter. The simulation of $H_S$ can be implemented in two different ways depending on the access model for the system Hamiltonian $H_S$. If $H_S$ is provided as a block-encoding then there exist optimal simulation techniques that add only a single extra ancilla qubit @low2019hamiltonian. If $H_S$ is provided as a sum of $k$-local Pauli strings or as a sparse matrix then product formula techniques can be used @childs2021theory at zero extra overhead in ancilla qubits. To simulate the time evolution of $alpha G$ we can use the following breakdown $e^(i alpha G t) = e^(i alpha U_("haar") D U_("haar")^dagger t) = U_("haar") e^(i alpha D t)U_("haar")^dagger$. This unitary can be implemented with a 2-design to approximate the $U_("haar")$, this is due to the fact that we only expand the channel to second order in $alpha$. Random Clifford circuits \cite{webb2015clifford} are sufficient for this purpose. To simulate $e^(i alpha D t)$ random $Z$ rotations and controlled-$Z$ rotations should be sufficient, as we only ever rely on the eigenvalues being pairwise independent in our analysis. In the worst case the number of random gates in the $Z$ basis that would need to be applied would scale with the dimension of the system, adding an overall factor of $dim_S$ to the preparation. In the product formula case we further remark that $H + alpha G$ could be simulated in total with a composite technique of using Trotter for $H$ and randomized compilation for $alpha G$ @hagan2023composite.

In classical Hamiltonian Monte Carlo it is well known that sharp gradients in the Hamiltonian require longer simulation time and more samples to converge. Our quantum routine has a much more subtle dependence on the structure of the Hamiltonian. As our single ancilla qubit only has one energy difference $gamma$, we have to tune this energy difference to allow for energy to be siphoned off from the system into the ancilla. This would present a conundrum, as knowing spectral gaps is as difficult or harder than preparing ground states of arbitrary quantum Hamiltonians, but we are able to prove that our routine is robust to complete ignorance of these differences. We show that this ignorance comes at an asymptotic cost in the amount of resources needed to prepare the thermal state. We numerically verify that knowledge of the eigenvalue differences can be used to speed up the total simulation time, as demonstrated in @fig_h_chain_error. We posit that this behavior serves as a crucial entry point for heuristics about Hamiltonian spectra into thermal state preparation algorithms. No prior thermal state preparation routines have had such an explicit demonstration of the utility of such knowledge. It was our hope to analytically quantify the speedups gained as a function of the relative entropy between a heuristic guess for the eigenvalue differences and the true spectra, but our numeric evidence will have to suffice until future work can clarify this dependence.

We would like to make a few remarks on potential improvements for the analysis of this channel. As we have demonstrated numerically, our guaranteed analytic values of $alpha$ and $t$ that lead to thermalization are drastically overestimated. We conjecture that this is due to our truncation of the weak-coupling expansion and in @fig_epsilon_scaling we demonstrate that taking $alpha prop 1 / t$ and $t prop 1 / sqrt(epsilon)$ drastically outperforms our analytically derived bounds of $alpha prop 1 / t^3$, by almost 4 orders of magnitude at $epsilon approx 0.005$. It is an open question of how to analyze this channel in the strong-coupling regime, and our numeric results suggest that such an analysis may indicate better performance of our protocol than a weak-coupling expansion can show. It is also an open question of whether dynamically chosen values of $alpha$ and $t$, such as having strong coupling and low time at the beginning and gradually decreasing $alpha$ and increasing $t$, can outperform static $alpha$ and $t$. We also suspect that the Markov relaxation theorem we used greatly overestimates the number of interactions needed. It remains to be seen if better Markov theory is needed or if the convergence time could be characterized based on the overlap of the initial state with the thermal state, which is a property that a few ground state preparation algorithms demonstrate. Another potential avenue for improving the analysis of this channel is whether different randomized interactions or even eigenvector heuristics can be beneficial. For example, in the harmonic oscillator if one has knowledge of the creation and annihilation operators $a^dagger$ and $a$, could one simply use the interaction $a^dagger tp (X + i Y) + a tp (X - i Y)$ instead of involving a randomized $G$ that relies on a Haar average? The last potential improvement is to extend our spectral gap computations using perturbation theory. We are only able to compute the spectral gap $tilde(lambda_star)(beta)$ in the limit of $beta -> oo$, but it should be possible to compute a perturbation on the order of $1 / beta$. This would give the simulation time needed to prepare low-temperature thermal states as opposed to zero-temperature states.

Lastly, we would like to speculate on possible applications of this routine to other quantum information processing tasks. The first question that arises is if these techniques could be used in the training of quantum Boltzmann machines, which are essentially thermal states. It is an open question if our thermalizing techniques could be used to either train models or to generate output samples from an already trained model. Through the process of demonstrating that this channel prepares the system in the thermal state we have calculated the output of our channel for both the system and the environment registers, and for much larger environments than single qubits. We can turn this protocol on it's head and ask how much information about the system are these ancilla qubits carrying away with them? Preliminary explorations suggest that given knowledge about eigenvalue gaps one can use transition statistics in the ancilla qubits to infer what the inverse temperature $beta$ is of the system, assuming the system is in a thermal state. Could this thermalizing channel instead be used to develop a Bayesian model to update beliefs about Hamiltonian spectra and system temperatures? This would represent an interaction agnostic model for performing quantum thermometry or spectroscopy, which to the best of our knowledge has not been developed yet.

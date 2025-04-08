// #import "macros.typ": *
#let ket(psi) = $lr(|#psi angle.r)$
#let bra(psi) = $lr(angle.l #psi|)$
#let ketbra(a, b) = $|#a angle.r angle.l #b|$
#let tp = $times.circle$
#let id = $bb(1)$


#import "@preview/ctheorems:1.1.3": *
#let lemma = thmplain("lemma", "Lemma", inset: (x: 0cm, top: 0cm))
#let proof = thmproof("proof", "Proof")
#let theorem = thmbox("theorem", "Theorem", stroke: 1pt, fill: rgb(0, 255, 0, 50))
#let definition = thmbox("definition", "Definition", stroke: 1pt, fill: rgb(0, 0, 255, 30))
#let todo = x => { text([TODO: #x], fill: red, weight: "bold") }

#show: thmrules.with(qed-symbol: $square$)
#import "conf.typ": *

#heading("Composite Simulations", level: 1, supplement: [Chapter]) <ch:composite_simulations>

The simulation of time-independent Hamiltonian dynamics is a fundamental primitive in quantum computing. To start, the computational problem of approximating the time dynamics of even $k$-local Hamiltonians (where $k$ is a small constant) is BQP-Complete. This means that any computational problem that can be solved efficiently on a quantum computer can be efficiently reduced to a simulation problem.

The simulation of quantum systems remains one of the most compelling applications for future digital quantum computers @whitfield2011simulation @jordan2012quantum @reiher2017elucidating @babbush2019quantum @su2021fault @o2021efficient.
As such, there are a plethora of algorithm options for compiling a unitary evolution operator $U(t) = e^{-i H t}$ to circuit gates @aharonov2003adiabatic @berry2007efficient @berry2015simulating @childs2019faster @low2019hamiltonian @low2019well @low2018hamiltonian @qdriftCampbell. Some of the simplest such algorithms are product formulas in which each term in a Hamiltonian $H = sum_i h_i H_i$ is implemented as $e^(i H_i t)$. A product formula is then a particular sequence
of these gates that approximates the overall operator $U(t)$. Two of the most well known product formula include Trotter-Suzuki Formulas @berry2007efficient @wiebe2010higher @childs2019faster @childs2021theory and the QDrift protocol in which terms are sampled randomly @qdriftCampbell @berry2020time. These two approaches are perhaps the most popular ancilla-free simulation methods yet discovered.

One of the main drawbacks of Trotter-Suzuki formulas is that each term in the Hamiltonian has to be included in the product formula regardless of the magnitude of the term. This leads
to a circuit with a depth that scales at least linearly with the number of terms in $H$, typically denoted $L$. QDrift avoids this by randomly choosing which
term to implement next in the product formula according to an importance sampling scheme in which higher weight terms have larger probabilities. The
downside to QDrift is that it has the same asymptotic scaling with $t / epsilon$ as a first-order Trotter formula, meaning it is outperformed at large
$t/ epsilon$ by even a second-order Trotter formula.

In this paper we present a framework for combining simulation channels in a way that allows one to flexibly interpolate the gate cost tradeoffs between the individual channels. The primary example we study is the composition of Trotter-Suzuki and QDrift channels. This is motivated in some part as an effort to extend
randomized compilers to include conditional probabilities and in some part to encapsulate progress in chemistry simulations of dropping small
weight terms or shuffling terms around different time steps @bucket_sim. This latter concept was first developed with the idea of "coalescing" terms into
"buckets" by Wecker et al. @bucket_sim and further explored by Poulin et al. @coalescing_con_wiebe. They showed that grouping terms of similar sizes together to be skipped during certain Trotter steps led to negligible increases in error and reduced gate counts by about a factor of 10. Similar improvements are also seen in the randomized setting of @kivlichan2019phase. In this work we extend on these ideas by placing a specific set of terms into a Trotter partition and the rest in a QDrift partition. This simple division can then be studied analytically and we are able to provide sufficient conditions on asymptotic improvements over completely Trotter or completely QDrift channels. Although we are not able to develop the idea of conditional samples in QDrift protocols, our
procedure can be viewed as a specific subset of what a generic Markovian QDrift would look like. We briefly mention these generalizations in
Section .

== Related Work

Recent approaches have sought to use the advantages of randomized compilation as a subset of an overall simulation, such as the hybridized scheme for interaction picture simulations @hybridized_interaction_pic. What separates these two works is that our approach offers a more flexible approach for generic time-independent simulation problems whereas the hybridized schemes are specifically tailored to taking advantage of the time dependence introduced by moving to an interaction picture. As such, the hybridized approach achieves asymptotic advantages when the size of the interaction picture term dominates the overall Hamiltonian. This typically occurs in instances in which the size of an operator is unbounded, which can occur in lattice field theory simulations or constrained systems. The way the hybridized scheme in @hybridized_interaction_pic works is via a "vertical" stacking of simulation channels, for example one channel to handle the Interaction Picture rotations and then other channels on top of this to simulate the time-dependence it generates on the remaining Hamiltonian terms. Our work instead remains in the Schrodinger time evolution picture and we perform a "horizontal" stacking of simulation techniques. By horizontal we mean for a given simulation time we split the Hamiltonian up into (potentially) disjoint partitions and simulate each partition for the full simulation time but with different techniques, such as Trotter or QDrift. These techniques allow us to achieve asymptotic improvements over either method for a loose set of assumptions.

There are two other simulation techniques that have been proposed recently that have a similar interpolation behavior between QDrift and Trotter channels. The first of these methods is the SparSto, or Stochastic Sparsification, technique by Ouyang, White, and Campbell @sparsto. The SparSto procedure randomly sparsifies the Hamiltonian and performs a randomly ordered first-order Trotter formula on the sampled Hamiltonian. They construct these probabilities such that the expected Hamiltonian is equal to the Hamiltonian being simulated. They then fix the expected number of oracle queries of the form $e^(i H_i t')$ and give diamond distance bounds on the resulting channel error. The claim for interpolation between Trotter and QDrift is that one can fix the expected number of gates to be 1 for each time step, in which case the sparsification mimics QDrift, whereas if no sparsification is performed then the channel is simply implementing Trotter. They show that this allows for one to have reduced simulation error up to an order of magnitude on numerically studied systems as compared to Trotter or QDrift. One downside to these techniques is that the number of gates applied is a random variable, so making gate cost comparisons is rather difficult especially considering that no tail bounds on high gate cost sampled channels are provided. In @sparsto they prefer to fix the expected gate cost and analyze the resulting diamond norm error. In contrast, our procedures directly implement both QDrift and Trotter channels and have a fixed, deterministic gate cost.

The second method of note with both QDrift and Trotter behavior is that of Jin and Li @jin2021partially. They develop an analysis of the variance of a unitary consisting of a first-order Trotter sequence followed by a QDrift channel. They focus on bounding the Mean Squared Error (MSE) of the resulting channel and use a simple partition of the Hamiltonian terms based on spectral norm. Their partitioning scheme places all terms below some cutoff into the first-order Trotter sequence and all terms above the cutoff into the QDrift channel. Their main results show an interpolation of the MSE between 0 when the partitioning matches a solely Trotter channel and matching upper bounds for QDrift when all terms are randomly sampled. This work goes beyond the results from Jin and Li by providing an analysis of the diamond distance between an ideal evolution and our implemented channel, which is more useful analytically than the MSE, as well as providing upper bounds on the number of gates needed in an implementation to meet this diamond distance. In addition our work remains independent of specific partitioning schemes as much as possible and instead places restrictions on which partitions achieve improvements. In the interest of practicality we do show methods for partitioning that can be useful in both the first-order and higher-order Trotter cases. Specifically for higher-order Trotter formulas we give a probabilistic partitioning scheme that is easily computable and matches gate cost upper bounds in the extreme limits as our probabilities saturate the QDrift and Trotter limits.

The rest of the paper is organized as follows. We first provide a brief summary of the main results in @sec:composite_main_results. After reviewing known results and notation in @sec:composite_prelim, we explore methods for creating Composite channels using First-Order Trotter Formulas with QDrift in @sec:composite_first_order as a warmup. This is broken down
into three parts in which we find the gate cost for an arbitrary partition, we then give a method for producing a good partitioning, and then we analyze conditions in which a Composite channel can beat either first-order Trotter or QDrift channels. In @sec:composite_higher_order we then extend this framework to more general higher-order Trotter Formulas. This section mirrors the organization of the first-order Trotter section,
namely we find the cost of an arbitrary partition, we give a method for producing a partition efficiently, and then we analyze when one could see
improvements over the constituent channels. Finally, in @sec:composite_discussion we discuss extensions to this model that allow a flexible interpolation between various types of product formulas that could be leveraged numerically.

== Main Results <sec:composite_main_results>

== Preliminaries <sec:composite_prelim>

In this section we will first introduce the necessary notation we will use and then state known results about Trotter-Suzuki formulas and QDrift channels. We work exclusively with time-independent Hamiltonians $H$ in a $2^n$ dimensional Hilbert space $cal(H)$. We also assume that $H$ consists of $L$ terms $H = sum_(i = 1)^L h_i H_i$ where $h_i$ represents the spectral norm of the term, $H_i$ is a Hermitian operator on $cal(H)$, and $norm(H_i) = 1$. Note without loss of generality we can always assume $h_i >= 0$, as we can always absorb the phase into the operator $H_i$ itself. We use $norm(M)$ to refer to the spectral norm, or the magnitude of the largest singular value of $M$. We use $lambda$ to refer to the sum of $h_i$, namely $lambda = sum_i h_i$. We will also use subscripts on $lambda$, such as $lambda_A$ to refer to sums of subsets of the terms of $H$. For example, if $H = 1 H_1 + 2 H_2 + 3 H_3$ and $G = 1 H_1 + 2 H_2$, then $lambda = 6$ and $lambda_G = 3$.

We use $U(t)$ to refer to the unitary operator $e^(i H t)$ and $cal(U)(t)$ to refer to the channel $rho |-> U(t) rho U(t)^dagger$. We will be particularly concerned with simulations of subsets of the terms of $H$, which we denote as follows. We typically work with a partition of $H$ into two matrices $H = A + B$, and we let $A = sum_i a_i A_i$ and $B = sum_j b_j B_j$, where we have simply relabeled the relevant $h_i$ and $H_i$ into $a$'s, $b$'s, $A$'s, and $B$'s. This allows us to define the exact unitary time evolution operators $U_A(t) = e^(i A t)$ and channels $cal(U)_A(t) = U_A(t) rho U_A(t)^dagger$, similarly defined for $B$. As we will be working with approximations to these channels, any operator or channel with a tilde represents an "implemented" channel, for example a first-order Trotter operator for $A$ would look like $tilde(U_A) (t) = e^(i a_1 A_1 t) dots e^(i a_L A_L t)$. We avoid using $cal(E)$ to represent an approximation or product formula as $cal(E)$ will be used for error channels.

Although much of the literature for Trotter-Suzuki formulas is written in terms of unitary operators $U = e^(i H t)$ acting on state vectors $ket(psi)$ for our purposes it will prove most natural to consider a product formula as a channel $cal(U) = e^(i H t) rho e^(-i H t)$ acting on a density matrix $rho$. After reviewing known results on unitary constructions of Trotter-Suzuki formulas we give a straightforward extension of these bounds to channels.

=== Product Formulas

We now show how to implement basic product formulas, namely Trotter-Suzuki or just Trotter formulas as well as QDrift, assuming access to arbitrary single qubit unitaries and controlled NOT gates. We will first show how to implement an arbitrary Pauli rotation $e^(i P t)$ for some Pauli string $P$ using no additional ancilla qubits. Then we will define the Trotter-Suzuki construction and give heuristic evidence for the first order scaling. We avoid giving a rigorous proof and instead refer the reader to the canonical paper by Childs et. al @childs2021theory. Lastly, we will present the construction of QDrift by Campbell @qdriftCampbell, providing a heuristic proof of correctness.

#definition("Trotter-Suzuki Formulae")[
    Given a Hamiltonian $H$, let $U_"TS"^((1))(t)$ denote the first-order Trotter-Suzuki time evolution operator
    $
        U_"TS"^((1))(t) := e^(i h_L H_L t) dots e^(i h_1 H_1 t) = product_(i = 1)^(L) e^(i h_i H_i t).
    $

    Note that the ordering of the factors in the product $product$ is defined to start from the rightmost operator and end at the leftmost: $product_(i = 1)^3 X_i = X_3 X_2 X_1$. Following this we can define the second-order Trotter-Suzuki time evolution operator as
    $
        U_"TS"^((2))(rho; t) &:= e^(i h_1 H_1 (t / 2)) dots e^(i h_L H_L (t / 2)) e^(i h_L H_L (t / 2)) dots e^(i h_1 H_1 (t / 2)) = product_(i = L)^(1) e^(i h_i H_i (t / 2)) product_(j = 1)^(L) e^(i h_j H_j (t / 2)).
    $

    This formula serves as the base case for the recursively defined higher-order formulas
    $ U_"TS"^((2k))(t) := U_"TS"^((2k-2))(u_k t)^2 dot U_"TS"^((2k-2))((1-4 u_k)t) dot U_"TS"^((2k-2)) (u_k t)^2, $
    where $u_k := 1 / (4-4^(1/(2k - 1)))$. In addition we define $Upsilon := 2 dot 5^(k-1)$ as the number of "stages" in the higher-order product formula. We can now introduce the time evolution channels as
    $ cal(U)_"TS"^((2k))(rho; t) := U_"TS"^((2k))(t) " " rho " " U_"TS"^((2k))(t)^dagger, $

    where for consistency we use the calligraphic $cal(U)^((2k))$ to represent the applied channels.
] <def:trotter_suzuki>
Despite their simplicity, Trotter-Suzuki formulas are fiendishly difficult to analyze. For decades the only error analysis that existed was worst-case analysis that often drastically overestimated the actual error. It was known that the first order expression depended on the commutator structure among the terms, but this was not generalized until 2021 in @childs2021theory, 25 years after Lloyd's original work @lloyd1996universal. We will follow @childs2021theory and denote the expression that captures this commutator scaling as $alpha_"comm"$ (in @childs2021theory $tilde(alpha)_"comm"$ is used) which is defined as
$
    alpha_"comm" (H, 2k) := sum_(gamma_i in {1, ..., L}) (product h_(gamma_i)) norm([H_(gamma_(2k + 1)), [H_(gamma_(2k)), [ ...,[H_(gamma_2), H_(gamma_1)] ... ]]])_oo .
$ <eq:alpha_comm_def>
We will also make heavy use of the restriction of $alpha_"comm"$ to subsets of the total Hamiltonian $H$. For example, if $H = A + B$ then we define the commutator structure over $A$ as
$
    alpha_"comm" (A, 2k) := sum_(gamma_i in {1, ..., L_A}) (product a_(gamma_i)) norm([A_(gamma_(2k + 1)), [A_(gamma_(2k)), [ ...,[A_(gamma_2),A_(gamma_1)] ... ]]])_oo .
$
This then allows us to decompose the total commutator structure into 3 pieces: commutators that contain only terms from $A$, commutators that contain only terms from $B$, and commutators that contain _at least_ one term from both $A$ and $B$
$
    alpha_"comm" (H, 2k) = alpha_"comm" (A, 2k) + alpha_"comm" (B, 2k) + alpha_"comm" ({A, B}, 2k).
$

This allows us to give the error associated with a a Trotter-Suzuki formula in the following theorem.
#theorem([Trotter-Suzuki Error @childs2021theory])[
    Let $U_"TS"^((2k))$ be the Trotter-Suzuki unitary as given in @def:trotter_suzuki for the Hamiltonian $H = sum_(i =1)^L h_i H_i$. Then the spectral norm of the difference between the $U_"TS"^((2k))(t/r)$ and the ideal evolution $U(t/r)$ is given by
    $
        norm(U(t\/r) - U_"TS"^((2k))(t\/r))_oo <= #math.cases(
        $(2 ( (Upsilon t )/ r)^(2k + 1) )/ (2k + 1)  alpha_"comm" (H, 2k) #h(9.5mm) " if " 2k > 1$,
        $t^2 / (2 r^2) " " alpha_"comm" (H,1) #h(20mm)" if " 2k = 1$
      )
    $
]
The complete proof of the above theorem is very nontrivial and beyond the scope of this thesis. See @childs2021theory for complete details, the proof of the higher order bounds can be found in Appendix E and the first order expression is found in Proposition 9 in Section V. Instead, we provide a heuristic proof for the first order error for completeness.

#proof("Heuristic First Order")[
    Compute a Taylor Series for the Trotter formula and the ideal evolution. First the ideal evolution:
    $
        U(t\/r) &= e^(i H t / r) = id + (i t) / r H + O((t\/r)^2).
    $
    Then the Trotter terms:
    $
        product_(i = 1)^L e^(i h_i H_i t') = product_(i = 1)^L (id + i h_i t' H_i + O(t'^2)) = id + i t' sum_(i = 1)^L h_i H_i + O(t'^2 ).
    $
    Pretty clear to see that in the difference $U(t/r) - U_"TS"^((1))(t/r)$ the zeroth and first order terms vanish, leaving only the second order.
]

=== Randomized Product Formulas
We now introduce QDrift @qdriftCampbell, one of the first randomized compilers for quantum simulation. The main idea that led to the development of QDrift is that instead of iterating through each term in the Hamiltonian to construct a product formula, or even a random ordering of terms as in @childs2019faster, each exponential is chosen randomly from the list of terms in $H$. Each term is selected with probability proportional to it's spectral weight, the probability of choosing $H_i$ is $h_i / (sum_j h_j) =: h_i / norm(h)$, and then simulated for a time $tau = norm(h) t$. This is the protocol for a single sample. As we will denote the portion of the Hamiltonian that we simulate with QDrift in later sections as $B$ we let $N_B$ denote the number of samples used.
#definition("QDrift Channel")[
    Let $N_B$ denote the number of samples, $norm(h) = sum_(i = 1)^L h_i$, and $tau := (norm(h) t)/ N_B$. The QDrift channel for a single sample is given as
    $
        cal(U)_"QD" (t; 1) := rho |-> sum_(i = 1)^L h_i / norm(h) e^(- i H_i norm(h) t) " " rho " " e^(+ i H_i norm(h) t),
    $
    and the QDrift channel for $N_B$ samples is
    $
        cal(U)_"QD" (t; N_B) := cal(U)_"QD" (t \/ N_B; 1)^(circle N_B).
    $
] <def:qdrift>
#h(1cm) Once we have the channel defined we can then state the main results of @qdriftCampbell.
#theorem("QDrift")[
    Given a Hamiltonian $H$, time $t$, and error bound $epsilon$, the ideal time evolution channel $cal(U)(t)$ can be approximated using $N_B = (4 t^2 norm(h)^2) / epsilon$ samples of a QDrift channel. This approximation is given by the diamond distance
    $
        norm(cal(U)(t) - cal(U)_"QD" (t; N_B))_diamond <= epsilon.
    $
    The number of operator exponentials is then given as
    $
        C_"QD" (H, t, epsilon) <= (4 t^2 norm(h)^2) / epsilon.
    $
]

== First Order Composite Channels <sec:composite_first_order>
We now turn towards combining the two product formulas given in @sec:composite_prelim in a Composite channel. We first will assume that the Hamiltonian has already been partitioned into two pieces $H = A + B$, where $A$ will be simulated with a first order Trotter formula and $B$ with QDrift. Given a fixed partitioning allows for us to compute the diamond distance error in the resulting channel, which then allows us to bound the number of operator exponentials needed to implement the channel. The resulting cost function will then be parametrized by the partitioning, which we can then use to determine an optimal partitioning algorithm. Finally, we give a specific instance in which a Composite channel can offer asymptotic improvements in query complexity over either a purely Trotter or QDrift channel.

=== Query Complexity <sec:composite_first_order_query_complexity>
To analyze the error of our Composite channel we need to first reduce the overall time evolution channel $rho |-> e^(-i H t) rho e^(+i H t)$ into the simpler pieces that we can analyze with our Trotter and QDrift results. Assuming a partitioning $H = A + B$, where $A$ consists of terms that we would like to simulate with Trotter and $B$ has the terms we would like to sample from with QDrift. We now introduce the "outer-loop" error $E_({A,B})$ induced by this partitioning, which is as follows
$
    E_({A,B})(t) := e^(-i H t) rho e^(+i H t) - e^(-i B t) e^(-i A t) rho e^(+i A t) e^(+i B t) .
$
We use the phrase "outer-loop" as this decomposition is done before any simulation channels are implemented.



== Higher Order Composite Channels <sec:composite_higher_order>

== Discussion <sec:composite_discussion>

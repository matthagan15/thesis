// #import "macros.typ": *
#let ket(psi) = $lr(|#psi angle.r)$
#let bra(psi) = $lr(angle.l #psi|)$
#let ketbra(a, b) = $|#a angle.r angle.l #b|$
#let tp = $times.circle$
#let id = $bb(1)$
#let dmd = $diamond.medium$

#set math.equation(number-align: bottom)

#import "@preview/ctheorems:1.1.3": *
#let lemma_og = thmbox("lemma", "Lemma", stroke: 1pt, bodyfmt: x => text(x, style: "italic"))
#let proof_og = thmproof("proof", "Proof", inset: (x: 0cm))
#let theorem_og = thmbox("theorem", "Theorem", stroke: 1pt, bodyfmt: x => text(x, style: "italic"))
#let definition_og = thmbox("definition", "Definition", stroke: 1pt, bodyfmt: x => text(x, style: "italic"))

// #import "@preview/theorion:0.3.3": *
// #import cosmos.simple
// #show: show-theorion

// #let theorem = thmbox("theorem", "Theorem", stroke: 1pt, fill: rgb(0, 255, 0, 50))
// #let definition = thmbox("definition", "Definition", stroke: 1pt, fill: rgb(0, 0, 255, 30))

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

We use $U(t)$ to refer to the unitary operator $e^(i H t)$ and $cal(U)(t)$ to refer to the channel $rho |-> U(t) rho U(t)^dagger$. We will be particularly concerned with simulations of subsets of the terms of $H$, which we denote as follows. We typically work with a partition of $H$ into two matrices $H = A + B$, and we let $A = sum_i a_i A_i$ and $B = sum_j b_j B_j$, where we have simply relabeled the relevant $h_i$ and $H_i$ into $a$'s, $b$'s, $A$'s, and $B$'s. This allows us to define the exact unitary time evolution operators $U_A(t) = e^(i A t)$ and channels $cal(U)_A(t) = U_A(t) rho U_A(t)^dagger$, similarly defined for $B$.

Although much of the literature for Trotter-Suzuki formulas is written in terms of unitary operators $U = e^(i H t)$ acting on state vectors $ket(psi)$ for our purposes it will prove most natural to consider a product formula as a channel $cal(U) = e^(i H t) rho e^(-i H t)$ acting on a density matrix $rho$. After reviewing known results on unitary constructions of Trotter-Suzuki formulas we give a straightforward extension of these bounds to channels.

=== Product Formulas

We now show how to implement basic product formulas, namely Trotter-Suzuki or just Trotter formulas as well as QDrift, assuming access to arbitrary single qubit unitaries and controlled NOT gates. We will first show how to implement an arbitrary Pauli rotation $e^(i P t)$ for some Pauli string $P$ using no additional ancilla qubits. Then we will define the Trotter-Suzuki construction and give heuristic evidence for the first order scaling. We avoid giving a rigorous proof and instead refer the reader to the canonical paper by Childs et. al @childs2021theory. Lastly, we will present the construction of QDrift by Campbell @qdriftCampbell, providing a heuristic proof of correctness.

#definition_og("Trotter-Suzuki Formulae")[
    Given a Hamiltonian $H$, let $S^((1))(t)$ denote the first-order Trotter-Suzuki time evolution operator and channel as
    $
        S^((1))(t) &:= e^(i h_L H_L t) dots e^(i h_1 H_1 t) = product_(i = 1)^(L) e^(i h_i H_i t), \
        cal(S)^((1))(rho; t) &:= S^((1))(t) dot rho dot S^((1))(t)^dagger.
    $ <eq:trotter_first_order>
    This first order formula serves as the base case for the recursively defined higher-order formulas
    #set math.equation(number-align: bottom)
    $
        S^((2))(t) &:= e^(i h_1 H_1 (t / 2)) dots e^(i h_L H_L (t / 2)) e^(i h_L H_L (t / 2)) dots e^(i h_1 H_1 (t / 2)) = S^((1))(t\/2) dot S^((1))(-t\/2)^dagger \
        S^((2k))(t) &:= S^((2k - 2))(u_k t) dot S^((2k - 2))(u_k t) dot S^((2k - 2)) ((1-4 u_k)t) dot S^((2k - 2))(u_k t) dot S^((2k - 2))(u_k t) \
        cal(S)^((2k))(rho;t) &:= S^((2k))(t) dot rho dot S^((2k))(t)^dagger,
    $ <eq:trotter_high_order>
    where $u_k := 1 / (4-4^(1/(2k - 1)))$. In addition we define $Upsilon := 2 dot 5^(k-1)$ as the number of "stages" in the higher-order product formula.
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
#theorem_og([Trotter-Suzuki Error @childs2021theory])[
    Let $S^((2k))$ be the Trotter-Suzuki unitary as given in @def:trotter_suzuki for the Hamiltonian $H = sum_(i =1)^L h_i H_i$. Then the spectral norm of the difference between the Trotter-Suzuki formulas $S^((1))(t\/r)$ and $S^((2k))(t\/r)$ and the ideal evolution $U(t\/r)$ is given by
    $
        norm(U(t\/r) - S^((1))(t\/r))_oo &<= t^2 / (2 r^2) " " alpha_"comm" (H,1), \
        norm(U(t\/r) - S^((2k))(t\/r))_oo &<= (( Upsilon t )^(2k + 1) ) / (r^(2k + 1)(k + 1\/2)) alpha_"comm" (H, 2k).
    $
] <thm_trotter_error>
The complete proof of the above theorem is very nontrivial and beyond the scope of this thesis. See @childs2021theory for complete details, the proof of the higher order bounds can be found in Appendix E and the first order expression is found in Proposition 9 in Section V. Instead, we provide a heuristic proof for the first order error for completeness.

#proof_og("Heuristic First Order")[
    Compute a Taylor Series for the Trotter formula and the ideal evolution. First the ideal evolution:
    $
        U(t\/r) &= e^(i H t / r) = id + (i t) / r H + O((t\/r)^2).
    $
    Then the Trotter terms:
    $
        product_(i = 1)^L e^(i h_i H_i t') = product_(i = 1)^L (id + i h_i t' H_i + O(t'^2)) = id + i t' sum_(i = 1)^L h_i H_i + O(t'^2 ).
    $
    Pretty clear to see that in the difference $U(t\/r) - S^((1))(t\/r)$ the zeroth and first order terms vanish, leaving only the second order.
]

=== Randomized Product Formulas
We now introduce QDrift @qdriftCampbell, one of the first randomized compilers for quantum simulation. The main idea that led to the development of QDrift is that instead of iterating through each term in the Hamiltonian to construct a product formula, or even a random ordering of terms as in @childs2019faster, each exponential is chosen randomly from the list of terms in $H$. Each term is selected with probability proportional to it's spectral weight, the probability of choosing $H_i$ is $h_i / (sum_j h_j) =: h_i / norm(h)$, and then simulated for a time $tau = norm(h) t$. This is the protocol for a single sample. As we will denote the portion of the Hamiltonian that we simulate with QDrift in later sections as $B$ we let $N_B$ denote the number of samples used.
#definition_og("QDrift Channel")[
    Let $N_B$ denote the number of samples, $norm(h) = sum_(i = 1)^L h_i$, and $tau := (norm(h) t)/ N_B$. The QDrift channel for a single sample is given as
    $
        cal(Q) (rho; t, 1) := sum_(i = 1)^L h_i / norm(h) e^(- i H_i norm(h) t) dot rho dot e^(+ i H_i norm(h) t),
    $
    and the QDrift channel for $N_B$ samples is
    $
        cal(Q) (t, N_B) := cal(Q) (t \/ N_B, 1)^(compose N_B).
    $
] <def:qdrift>
#h(1cm) Once we have the channel defined we can then state the main results of @qdriftCampbell.
#theorem_og("QDrift Cost", supplement: "QDrift Cost")[
    Given a Hamiltonian $H$, time $t$, and error bound $epsilon$, the ideal time evolution channel $cal(U)(t)$ can be approximated using $N_B = (4 t^2 norm(h)^2) / epsilon$ samples of a QDrift channel. This approximation is given by the diamond distance
    $
        norm(cal(U)(t) - cal(Q) (t, N_B))_(dmd) <= epsilon.
    $
    The number of operator exponentials is then given as
    $
        C_"QD" (H, t, epsilon) <= (4 t^2 norm(h)^2) / epsilon.
    $
] <thm:qdrift_cost>

#proof_og([@thm:qdrift_cost])[
    This one is also a Taylor's Series
    $
        cal(U)_"QD" (t)
    $
]

== First Order Composite Channels <sec:composite_first_order>
We now turn towards combining the two product formulas given in @sec:composite_prelim in a Composite channel. We first will assume that the Hamiltonian has already been partitioned into two pieces $H = A + B$, where $A$ will be simulated with a first order Trotter formula and $B$ with QDrift. Given a fixed partitioning allows for us to compute the diamond distance error in the resulting channel, which then allows us to bound the number of operator exponentials needed to implement the channel. The resulting cost function will then be parametrized by the partitioning, which we can then use to determine an optimal partitioning algorithm. Finally, we give a specific instance in which a Composite channel can offer asymptotic improvements in query complexity over either a purely Trotter or QDrift channel.

=== Query Complexity <sec:composite_first_order_query_complexity>
To analyze the error of our Composite channel we need to first reduce the overall time evolution channel $rho |-> e^(-i H t) rho e^(+i H t)$ into the simpler pieces that we can analyze with our Trotter and QDrift results. Assuming a partitioning $H = A + B$, where $A$ consists of terms that we would like to simulate with Trotter and $B$ has the terms we would like to sample from with QDrift. We now introduce the "outer-loop" error $E_({A,B})$ induced by this partitioning, which is as follows
$
    E_({A,B})(t) := e^(-i H t) rho e^(+i H t) - e^(-i B t) e^(-i A t) rho e^(+i A t) e^(+i B t) .
$
We use the phrase "outer-loop" as this decomposition is done before any simulation channels are implemented. This gives the first order Composite channel $cal(C)$ as
$
    cal(C)^((1))(t) &:= cal(Q)_B (t, N_B) compose cal(S)_A^((1))(t).
$
We will first bound the error of this approximation to the ideal evolution. This error bound will then allow us to bound the number of exponentials needed to approximate the ideal dynamics.
#lemma_og("First-Order Composite Error")[
    Given a Hamiltonian $H$ partitioned into a first order Trotter term $A$ and QDrift term $B$ such that $H = A + B$, the first order Composite Channel $cal(C)^((1))$ has an error of at most
    $
        norm(cal(U)(t) - cal(C)^((1))(t) )_(dmd) <= t^2 (sum_(i,j) a_i a_j norm([A_i, A_j]) + sum_(i,j) norm([A_i, B_j ]) + (4 norm(h_B)^2) / N_B).
    $
] <lem:composite_error_first_order>
#proof_og()[
    #set math.equation(number-align: bottom)
    We will first use an outer-loop decomposition to get the error associated by our partitioning, note we temporarily supress arguments of $(t)$ for clarity,
    $
        norm(cal(U) - cal(C)^((1)))_dmd &= norm(cal(U) - cal(U)_B compose cal(U)_A  + cal(U)_B  compose cal(U)_A  - cal(C)^((1)) ) \
        &<= norm(cal(U) - cal(U)_B compose cal(U)_A)_dmd + norm(cal(U)_B compose cal(U)_A - cal(C)^((1)))_dmd.
    $<tmp_composite_0>
    We then can bound the leftmost term using the error decomposition
    $
        norm(cal(U) - cal(U)_B compose cal(U)_A)_dmd <= 2 norm(U - U_B dot U_A)_oo <= t^2 sum_(i,j) a_i b_j norm([A_i, B_j])_oo.
    $ <tmp:composite_1>
    And the rightmost term can be bounded using the subadditivity of the diamond distance

    $
        norm(cal(U)_B compose cal(U)_A - cal(C)^((1)) )_dmd &= norm(cal(Q)_B compose cal(S)^((1))_A - cal(U)_B compose cal(U)_A)_dmd \
        &<= norm(cal(U)_B - cal(Q)_B)_dmd + norm(cal(U)_A - cal(S)_A^((1)))_dmd \
        &<= (4 norm(h_B) t^2) / N_B + t^2 sum_(i,j) a_i a_j norm([A_i, A_j])_oo.
    $ <tmp:composite_2>
    Substituting @tmp:composite_1 and @tmp:composite_2 into @tmp_composite_0 yields the statement.
]
#h(5mm) Now that we have an upper bound on first order error for an arbitrary $t$ we can leverage this into a bound on the number of operator exponentials to reach arbitrary error $epsilon$ using standard time-slicing arguments. By letting $t -> t\/r$ and then repeating our Composite channel $r$ times we can control the accumulated error from each step. One of the beautiful features of product formulas is that this time-slicing leads to an overall reduction in the error, or in other words
$
    lim_(r -> oo) (e^(i A t / r) e^(i B t / r))^r = e^(i (A + B) t).
$
The following theorem utilizes a quantitative variant of the above, along with the error bounds we just proved, to provide the first order Composite cost Theorem.
#theorem_og("First-Order Composite Cost")[
    Given a time $t$, error bound $epsilon$, and a partitioned Hamiltonian $H = A + B$, the first order Composite Channel $cal(C)^((1))$ approximates the ideal time evolution operator $norm(cal(U)(t) - cal(C)^((1))(t\/r)^(compose r))_dmd <= epsilon$ using no more than
    $
        C_"Comp"^((1)) = (L_A + N_B) ceil(t^2/epsilon (sum_(i,j) a_i a_j norm([A_i, A_j]) + sum_(i,j) norm([A_i, B_j ]) + (4 norm(h_B)^2) / N_B))
    $
    operator exponential queries.
] <thm_composite_first_order_cost>
#proof_og()[
    We first will use the fact that since $H$ commutes with itself the time evolution operator can be decomposed into $r$ steps as $cal(U)(t) = cal(U)(t\/r)^(compose r)$. Then we can use the sub-additivity of the diamond norm with respect to channel composition to get the bound
    $
        norm(cal(U)(t) - cal(C)^((1))(t\/r)^(compose r))_dmd &= norm(cal(U)(t\/r)^(compose r) - cal(C)^((1))(t\/r)^(compose r))_dmd <= r norm(cal(U)(t\/r) - cal(C)^((1))(t\/r))_dmd .
    $
    Now we utilize @lem:composite_error_first_order to upper bound the single step error
    $
        norm(cal(U)(t) - cal(C)^((1))(t\/r)^(compose r))_dmd <= r (t^2 / r^2)(sum_(i,j) a_i a_j norm([A_i, A_j]) + sum_(i,j) norm([A_i, B_j ]) + (4 norm(h_B)^2) / N_B)."     "
    $
    In order for the above to be upper bounded by $epsilon$ we require
    $
        r >= (t^2 / epsilon)(sum_(i,j) a_i a_j norm([A_i, A_j]) + sum_(i,j) norm([A_i, B_j ]) + (4 norm(h_B)^2) / N_B),
    $
    and since increasing $r$ only increases the number of operator exponentials used we simply set $r$ to be the ceiling of the right hand side. This then yields the theorem as we have $r$ applications of $cal(C)^((1))$ and each application uses $L_A$ operator exponentials for the Trotter channel and $N_B$ samples of the QDrift channel.
]

=== First-Order Parameter Settings <sec_first_order_parameter_settings>

Our next task will be to determine "parameter settings" that optimize this gate cost, namely the partitioning $A + B$ and setting the number of QDrift samples $N_B$. To do this it will prove useful to have a continuous, non-integer variant of the gate cost expression which we define as
$
    tilde(C)_"Comp"^((1)) := (L_A + N_B) t^2 / epsilon (sum_(i,j) a_i a_j norm([A_i, A_j]) + sum_(i,j) norm([A_i, B_j ]) + (4 norm(h_B)^2) / N_B).
$
This is the same as $C_"Comp"^((1))$ but without the ceiling operation $ceil(dot)$. Although a user could use any value of $N_B$ they want, such as always setting $N_B = 1$, we provide the following setting for $N_B$ that is optimal with respect to the continuous variant of the gate cost.
#lemma_og()[
    Let $tilde(C)_"Comp"^((1))$ denote the continuous relaxtion to the cost of a first-order Composite channel with a given partitioning $H = A + B$. The optimal assignment of the number of QDrift samples $N_B$ with respect to $tilde(C)_"Comp"^((1))$ is given by
    $
        N_B = (2 norm(h_B) sqrt(L_A)) / sqrt(  (sum_(i,j) a_i a_j norm([A_i, A_j]) + sum_(i,j) a_i b_j norm([A_i, B_j]) )).
    $ <eq_composite_first_order_nb>
    This assignment is not valid if both $[A_i, A_j] = 0$ for all $A_i, A_j$ and $[A_i, B_j] = 0$ for all $A_i$ and $B_j$.
] <lem_composite_first_order_optimal_nb>
#proof_og()[
    We first compute the derivative of $tilde(C)_"Comp"^((1))$ with respect to $N_B$ as
    $
        (diff tilde(C)_"Comp") / (diff N_B) = t^2 / epsilon [sum_(i,j) a_i a_j norm([A_i, A_j]) + sum_(i,j) a_i b_j norm([A_i, B_j]) - (4 norm(h_B)^2 L_A) / N_B^2].
    $
    Setting this equal to 0 and solving for $N_B$ yields @eq_composite_first_order_nb. We then compute the second derivative as
    $
        (diff^2 tilde(C)_"Comp"^((1))) / (diff N_B^2) = (4 t^2 norm(h_B)^2 L_A) / (epsilon N_B^3),
    $
    which is always positive and therefore indicates that the optima found is a minimal cost with respect to $N_B$.
]

#h(5mm) There is still one remaining problem with the first-order Composite channel that we must address before we can compare it to existing product formulas: partitioning. Up until now we have assumed that a partitioning was given, but this is not a realistic assumption to make. There are many heuristics that one could, and most likely should, take advantage of when implementing an actual channel. For example, in a chemistry simulation one can put the nuclei-electron interactions, which are typically stronger than the electron-electron interactions, into Trotter and use a spectral norm cutoff to determine the remainder. On could construct a Hubbard like model on a grid but with long-range interactions, treating the "tunneling" kinetic term and nearest neighbor interactions with Trotter and then sampling the long-range interactions with QDrift. These kinds of heuristics will most likely be important in reducing simulation costs but will ultimately be application dependent.

We would like to provide a general purpose algorithm, one that works for any Hamiltonian and matches our intuition that the above heuristics capture. The algorithm we propose is a gradient based scheme that is based on a weighting of the Hamiltonian terms $H_i$. This weighting is based on the following trick, for every term $H_i$ we introduce a parameter $w_i$ that accounts for the weight of $H_i$ in the Trotter partition. Then our Hamiltonian can be written as
$
    H = sum_i h_i H_i = sum_i (w_i h_i H_i + (1- w_i) h_i H_I) = sum_i w_i h_i H_i + sum_i (1 - w_i) h_i H_i.
$
Then our partitioning is just $A = sum_i w_i h_i H_i$ and $B = sum_i (1 - w_i) h_i H_i$. We also would like to keep our interpretation of $w_i h_i$ and $(1 - w_i) h_i$ as spectral norms of the respective term in $A$ and $B$, so we restrict $w_i in [0,1]$. In this sense the weights could be thought of as probabilities, but we will not make use of any expectations or other probabilistic notions here so they should just be thought of as weight parameters. We will make use of a probabilistic variant of this scheme in @sec:composite_higher_order.

Once an initial weighting of each term is chosen, say $w_i = 1$ for the term with the largest spectral norm $h_"max"$ and $w_i = 0$ for every other term or maybe $w_i = 1/2$ for all $i$, we propose a greedy algorithm for computing a new set of weights $w_i '$. This greedy algorithm is based on the following gradient calculation of the weighted first order Composite channel cost
$
    (diff tilde(C)_"Comp"^((1))) / (diff w_m) = (L_A + N_B) t^2 / epsilon (h_m sum_j h_j norm([H_j, H_m]) - (8 h_m sum_i (1 - w_i) h_i) / N_B ) .
$
This gradient could be used in gradient descent, once the derivatives for all $w_m$ are computed and put into a vector one could update the parameters with $w_i ' = w_i - eta gradient_w tilde(C)_"Comp"^((1))$, where $eta$ is some learning rate.

Although this gradient descent based algorithm would be relatively easy to compute and implement in practice, it does not take advantage of the fact that our gradients are not only analytic, but linear with respect to a _single_ parameter $w_m$. This means that if we only update a single weight parameter at a time we can find an optimal assignment of $w_m$ with respect to the partial derivative given above. Setting the derivative equal to 0 and solving for $w_m$ yields
$
    (diff tilde(C)_"Comp"^((1)) ) / (diff w_m) = 0 ==> w_m ' = 1 - sum_(i != m) h_i / h_m ( norm([H_i, H_m]) / 8 - (1 - w_i)).
$ <eq_composite_first_order_partition>
Once $w_m '$ is determined then we can move on to $w_(m + 1) '$ until all parameters have been updated, at which point we can repeat until the parameters do not change. We unfortunately cannot provide more analysis, such as how many iterations this process will take, due to the coupling of the parameters.

Although our expression does not give exact partitionings that can be calculated in one go, our greedy algorithm based on @eq_composite_first_order_partition does capture two key pieces of intuition that we suspect good partitions will have. The first is that if a term $H_m$ commutes with every other term ($forall i != m " " [H_i, H_m] = 0$) then we have that $w_m ' = 1 + sum_(i != m) (h_i (1 - w_i)) / h_m >= 1,$ so we then restrict the update to $w_m ' = 1$. This implies that $H_m$ is completely placed into the Trotter partition, as we would expect.
The other piece of intuition is that smaller terms are pushed more towards the QDrift side of the partitioning. This can be seen from @eq_composite_first_order_partition while considering the limit as $h_m -> 0$. If we assume that $norm([H_i, H_m]) >= (1-w_i)$ on average, then the expression becomes $w_m -> -oo$ in this limit, which we stop at 0. In total, this indicates that large terms that do not commute with small terms with most of their weight in the Trotter partition tend to push the weight of small terms more towards QDrift.

=== Comparison with Trotter and QDrift <sec_composite_first_order_comparison>

Now that we have analyzed the cost and given a partitioning scheme for the first order Composite channel, we would like to know under what conditions this Composite channel can lead to comparable errors with lower gate cost. Instead of aiming to show that a Composite channel will outperform either first-order
Trotter or QDrift for arbitrary Hamiltonians, we instead illustrate a concrete setting in which we achieve guaranteed asymptotic improvements. In
later sections we are able to show more generic conditions in which asymptotic improvements can be obtained for higher-order formulas.

To start, let $H$ be a Hamiltonian that has a partitioning into $A$ and $B$ such that the following conditions hold.
+ The number of non-zero commutators between terms in $A$ scales with the square root of $L_A$. Mathematically, $ |{(i,j): norm([A_i, A_j]) != 0}| =: N_"nz"^2 in o(L_A). $
+ The strength of the $B$ terms, $norm(h_B) = sum_i b_i$, is asymtotically less than the maximum commutator norm divided by the number of terms in $A$ $ norm(h_B) <= (a_"max" N_"nz"^2) / L_A,  $ <tmp_composite_3> where $a_"max" = max_i a_i$.
+ The number of terms in the $A$ partition is vanishingly small compared to the total number of terms: $L_A in o(L)$.
Next, we can use the optimal $N_B$ value from @lem_composite_first_order_optimal_nb and @tmp_composite_3 to show that
$
    N_B^(-1) in O(1 / norm(h_B) sqrt(N_"nz"^2 a_"max"^2 + L_A a_"max" norm(h_B))) = O((N_"nz" a_"max") / norm(h_B) ).
$
Similarly we have $N_B in O( (norm(h_B) sqrt(L_A) ) / (a_"max" N_"nz"))$. Thus, @thm_composite_first_order_cost gives us the asymptotic number of operator exponentials as
$
    C_"Comp"^((1)) &in O(t^2 / epsilon (L_A + (norm(h_B) sqrt(L_A)) / (a_"max" N_"nz") )(a_"max"^2 N_"nz"^2 + L_A a_"max" norm(h_B) + (norm(h_B) N_"nz" a_"max") / sqrt(L_A) ) ) \
    &in O((t^2 L_A) / epsilon (a_"max"^2 N_"nz"^2)) \
    &in o(t^2 / epsilon L_A^2 a_"max"^2).
$
The lowest order Trotter formula for this simulation has the following asymptotic operator exponential cost, as given in
$
    C_"Trot"^((1)) in O(t^2 / epsilon (L N_"nz"^2 a_"max"^2)) in o(t^2 / epsilon L L_A a_"max"^2) subset.eq omega(C_"Comp"^((1))).
$
For QDrift we can use @thm:qdrift_cost to compute
$
    C_"QD" &in O(t^2 / epsilon (L_A a_"max" + norm(h_B))^2) \
    &in O((t^2 L_A^2 a_"max"^2) / epsilon (1 + (N_"nz" / L_A)^2)) \
    &in O((t^2 L_A^2 a_"max"^2) / epsilon) \
    &in omega(C_"Comp"^((1))).
$
This shows that there exist circumstances where the cost of the Composite channel scales better than either of the two methods that compose it, i.e. $C_"Comp"^((1)) in o(min(C_"Trot"^((1)), C_"QD"))$.

== Higher Order Composite Channels <sec:composite_higher_order>

We now move on from first-order Trotter formulas to arbitrary higher-order Trotter formulas. To analyze this case there are a few distinct differences with the first-order channels. The first is that we now have a choice for what order formula we would like to use for the outer-loop decomposition. Previously, for the first-order decomposition we used $cal(U)_B compose cal(U)_A$, but it will prove useful in our analysis to match the outer-loop order with the inner-loop Trotter order. For example, a second order outer-loop decomposition would look like $cal(U)_A (t\/2) compose cal(U)_B (t) compose cal(U)_A (t\/2)$, where we combined the two innermost $cal(U)_B (t\/2)$ for compactness. The next difference is that the time scaling between QDrift, Trotter, and the outer-loop errors could all be of different orders in $t/r$ which leads to a non-analytically solvable polynomial in $r$. The last issue that we address is that the commutator structure is no longer quadratic with respect to the Hamiltonian spectral norms, so we cannot follow the term weighting partitioning scheme from the first-order case. We will follow the same organizational structure as the first-order case and first set up our definitions and bound the diamond distance error, then compute the number of $e^(i H_i t)$ queries, followed by developing a partitioning scheme, and finally discuss the cost comparisons between our Composite channel and its constituents.

=== Query Complexity
In order to determine the number of queries needed for a Composite channel to approximate $cal(U)$ we first need to bound the diamond distance error for a single iteration. We will then use time-slicing arguments similar to the proof of @thm_composite_first_order_cost to compute the number of operator exponentials required for an accurate approximation. First, we need to give a rigorous definition of the higher order Composite channel.
#definition_og([Higher Order Composite Channel])[
    Given a Hamiltonian $H$ partitioned into two terms $A$ and $B$, let $cal(C)^((2k, 2l))$ denote the associated Composite channel that utilizes a $2k^"th"$ order inner-loop for the Trotter-Suzuki partition $cal(S)_A^((2k))$ and has a $2l^"th"$ order outer-loop. The outer-loop construction for the Composite channel can be constructed recursively from the base case for $l = 1$, which is given by
    $
        cal(C)^((2k, 2))(t) := cal(Q)_B (t\/2) compose cal(S)_A^((2 k)) (-t\/2)^dagger compose cal(S)_A^((2 k)) (t\/2) compose cal(Q)_B (t\/2),
    $
    and the higher-order outer-loop Composite channels are recursively defined as
    $
        cal(C)^((2k, 2l))(t) := cal(C)^((2k, 2l - 2)) (u_l t)^(compose 2) compose cal(C)^((2k, 2l - 2))((1-4 u_l)t) compose cal(C)^((2k, 2l - 2)) (u_l t)^(compose 2),
    $ <eq_composite_higher_order_def>
    where $u_l$ and the number of stages $Upsilon$ are the same as in @def:trotter_suzuki. We will typically ignore the distinction between inner and outer loops and use $cal(C)^((2k)) = cal(C)^((2k, 2k))$.
]

#lemma_og("Higher-Order Composite Error")[
    The diamond distance of a single higher-order Composite channel to the ideal time evolution channel is upper bounded as
    $
        norm(cal(U)(t) - cal(C)^((2k))(t))_dmd <= 2 (Upsilon t)^(2k + 1) / (k + 1 \/ 2) (alpha_"C" ({A, B}, 2k) + Upsilon alpha_"C" (A, 2k)) + Upsilon (4 norm(h_B)^2 t^2) / N_B.
    $
    We will use the following definitions for brevity
    $
        P(t) &:= t^(2k + 1) (2 Upsilon^(2k + 1)) / (k + 1\/2) (alpha_"comm" ({A, B}, 2k) + alpha_"comm" (A, 2k)) \
        Q(t) &:= t^2 (4 Upsilon norm(h_B)^2) / N_B.
    $ <eq_composite_p_n_q_def>
] <lem_composite_higher_order_error>
#proof_og([of @lem_composite_higher_order_error])[
    $
        norm(cal(U)(t) - cal(C)^((2k))(t))_dmd &<= norm(cal(U)(t) - cal(S)^((2k))({A, B}, t))_dmd + norm(cal(S)^((2k))({A,B}, t) - cal(C)^((2k))(t))_dmd \
        &<= 2 norm(e^(i H t) - S^((2k))({A, B}, t)) + norm(cal(S)^((2k))({A,B}, t) - cal(C)^((2k))(t))_dmd.
    $
    We can use @thm_trotter_error to bound the outer-loop error on the left as
    $
        norm(e^(i H t) -  S^((2k))({A, B}, t)) <= (Upsilon t)^(2 k + 1) / (k + 1\/ 2) alpha_"comm" ({A, B}, 2k).
    $
    We then use an inductive proof to argue that the inner-loop errors can be bounded as
    $
        norm(cal(S)^((2l))({A,B}, t) - cal(C)^((2k, 2l))(t))_dmd <= Upsilon (norm(cal(U)_A (t) - cal(S)_A^((2k))(t))_dmd + norm(cal(U)_B (t) - cal(Q)_B (t))_dmd),
    $ <tmp_composite_4>
    where the induction is over the outer-loop indexing of $2l$.
    - *Base Case ($2l = 1$):* $
        &norm(cal(C)^((2k,2))(t) - cal(S)^((2))({A, B})(t))_dmd \
        =& norm( cal(Q)_B (t\/2) compose cal(S)_A^((2k))(-t\/2)^dagger compose cal(S)_A^((2k))(t\/2) compose cal(Q)_B (t\/2) - cal(U)_B (t\/2) compose cal(U)_A (-t\/2)^dagger compose cal(U)_A (t\/2) compose cal(U)_B ( t\/ 2) )_dmd \
        <=& 2 norm(cal(U)_A (t\/2) - S_A^((2)) (t\/2))_dmd + 2 norm(cal(U)_B (t\/2) - cal(Q)_B (t\/2))_dmd.
    $ Since $Upsilon_1 = 2$ this matches the induction hypothesis.
    - *Inductive Step: * In this scenario we assume that the hypothesis in @tmp_composite_4 holds for $2l - 2$ and we would like to show it holds for $2l$. We do so via the recursive structure given in @eq_composite_higher_order_def and @def:trotter_suzuki, which allows us to express the hypothesis as
    $
        &norm(cal(C)^((2k, 2l))(t) - cal(S)^((2l))({A,B}, t))_dmd \
      =& norm( cal(C)^((2k, 2l - 2)) (u_l t)^(compose 2) compose cal(C)^((2k, 2l-2))((1-4 u_l) t) compose cal(C)^((2k, 2l - 2)) (u_l t)^(compose 2) \
      &  - cal(S)^((2l - 2))({A,B}, u_l t)^(compose 2) compose cal(S)^((2l -2))({A, B}, (1 - 4 u_l) t) compose cal(S)^((2l - 2))({A,B}, u_l t)^(compose 2))  #place($diamond.small$, dy: +1mm, dx: -0mm) \
      <=& 4 norm(cal(C)^((2k, 2l -2))(u_l t) - cal(S)^((2l - 2))({A,B}, u_l t))_dmd \
      &" " +norm(cal(C)^((2k, 2l -2))((1 - 4 u_l) t) - cal(S)^((2l - 2))({A,B}, (1 - 4 u_l) t))_dmd \
      <=& 4 Upsilon_(l - 1) (norm(cal(U)_A (t) - cal(S)_A^((2k))(t))_dmd + norm(cal(U)_B (t) - cal(Q)_B (t))_dmd) \
      &+ Upsilon_(l - 1) (norm(cal(U)_A (t) - cal(S)_A^((2k))(t))_dmd + norm(cal(U)_B (t) - cal(Q)_B (t))_dmd) \
      &= 5 Upsilon_(l - 1) (norm(cal(U)_A (t) - cal(S)_A^((2k))(t))_dmd + norm(cal(U)_B (t) - cal(Q)_B (t))_dmd) \
      =& Upsilon_l (norm(cal(U)_A (t) - cal(S)_A^((2k))(t))_dmd + norm(cal(U)_B (t) - cal(Q)_B (t))_dmd).
    $
    Therefore the inductive step holds.

    One point of emphasis we would like to make is that we are explicitly not utilizing the smaller times that come with the recursive outer-loop step. This is simply due to the difficulty of bookkeeping for each different time step used and it will be sufficient to use the upper bound of $t$. We can now continue with our proof of the lemma by substituting in the known Trotter-Suzuki and QDrift error terms, leading us to
    $
        norm(cal(U)(t) - cal(C)^((2k))(t))_dmd &<= 2 (Upsilon t)^(2k + 1) / (k + 1 \/ 2) alpha_"comm" ({A, B}, 2k) + Upsilon (norm(cal(U)_A (t) - cal(S)_A^((2k))(t))_dmd + norm(cal(U)_B - cal(Q)_B (t))_dmd ) \
        &<= 2 (Upsilon t)^(2k + 1) / (k + 1 \/ 2) (alpha_"comm" ({A, B}, 2k) + Upsilon alpha_"comm" (A, 2k)) + Upsilon (4 norm(h_B)^2 t^2) / N_B .
    $
    Note that the extra factor of $Upsilon$ in front of $alpha_"comm" (A)$ comes from the fact that we have $Upsilon$ copies of the $A$ simulation channel as opposed to just one outer-loop decomposition.
    // It will prove advantageous for us to have simple expressions for these two errors, so we define the following
    // $
    //     P(t) &:= t^(2k + 1) (2 Upsilon^(2k + 1)) / (k + 1\/2) (alpha_"comm" ({A, B}, 2k) + alpha_"comm" (A, 2k)) \
    //     Q(t) &:= t^2 (4 Upsilon norm(h_B)^2) / N_B,
    // $
    // where "$P$" can be thought of as "product formula" and $Q$ for QDrift.
]

#h(5mm) Now that we have bounded the diamond distance error for a single time step we can proceed with our time-slicing arguments to produce a controllable error bound. This will lead us to our final expression for the query cost of a higher order composite method.
#theorem_og("Higher-Order Composite Cost")[
    Given a time $t$, error bound $epsilon$, and a partitioned Hamiltonian $H = A + B$ the $2k^"th"$ order Composite channel $cal(C)^((2k))$ utilizes at most
    $
        C_"Comp"^((2k)) <= Upsilon (Upsilon L_A + N_B) ceil( (4^(1/(2k)) (Upsilon t)^(1 + 1/2k)) / ((2k + 1) epsilon)^(1/(2k)) (alpha_"C" ({A,B}) + Upsilon alpha_"C" (A)) + (4 Upsilon norm(h_B)^2 t^2) / (N_B epsilon) ),
    $
    gates,where the $alpha_"C"$ are both of order $2k$, to meet the error budget given by
    $
        norm(cal(U)(t) - cal(C)^((2k))(t\/r)^(compose r))_dmd <= epsilon.
    $
    Using the upper bounds provided for Trotter-Suzuki and QDrift evolution channels and defining
    $
        q_B := (alpha_"Comm" (B, 2k) ) / (alpha_"Comm" (H, 2k))
    $
    to capture the amount of "commutator structure" of $H$ that is contained in $B$, we can rewrite this upper bound as
    $
        C_"Comp"^((2k)) <= Upsilon (Upsilon L_A + N_B) ceil((C_"Trot"^((2k))(H, t, epsilon) )/ Upsilon^(1 - 1\/2k) (1 - q_B)^(1\/2k) / ( L) + C_"QD" (H, t, epsilon) Upsilon / N_B (norm(h_B) / norm(h))^2).
    $
]
#proof_og()[
    We start off by utilizing standard time-slicing arguments to bound our single step distance
    $
        norm(cal(U)(t) - cal(C)^((2k))(t\/r)^(compose r))_dmd = norm(cal(U)(t\/r)^(compose r) - cal(C)^((2k))(t\/r)^(compose r))_dmd <= r norm(cal(U)(t\/r) - cal(C)^((2k))(t\/r))_dmd.
    $
    Using our results in @lem_composite_higher_order_error we can then bound the single time step error as
    $
        norm(cal(U)(t\/r) - cal(C)^((2k))(t\/r))_dmd <= P(t) / r^(2k + 1) + Q(t) / r^2.
    $

]

== Numerics

== Discussion <sec:composite_discussion>
Testing theorion

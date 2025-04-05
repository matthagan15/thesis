// #import "macros.typ": *
#let ket(psi) = $lr(|#psi angle.r)$
#let bra(psi) = $lr(angle.l #psi|)$
#let ketbra(a, b) = $|#a angle.r angle.l #b|$

#let tp = $times.circle$
#let id = $bb(1)$


#import "@preview/ctheorems:1.1.3": *
#let lemma = thmplain("lemma", "Lemma", inset: (x: 0cm, top: 0cm))
#let proof = thmproof("proof", "Proof")
#let definition = thmplain("definition", "Definition")
#show: thmrules.with(qed-symbol: $square$)
#import "conf.typ": *

// = Composite Simulations <ch:composite_simulations>
#heading("Composite Simulations", level: 1, supplement: "Chapter") <ch:composite_simulations>

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

=== Main Results <sec:composite_main_results>

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

== First Order Composite Channels <sec:composite_first_order>

== Higher Order Composite Channels <sec:composite_higher_order>

== Discussion <sec:composite_discussion>

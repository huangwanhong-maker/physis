# 0020. Chains on a Classified Set, and the Probability of a Record of Events

- **Status:** Accepted, provisionally. Cells change independently between pieces, and kernels are dense over the cells they join.
- **Date:** 2026-11-04
- **Decided by:**
  - **The request.** The user asked for the stages E3 and E4 of the design of the infrastructure's Event Spacetime ([ADR 0067](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/development/decisions/0067-event-spacetime.md)): "Then we can begin E3 and E4".
  - **What Claude did.** It wrote these two headers and this record.

## 1. Context

**The design.** The design of Event Spacetime ([development/ideas/2026-11-04-event-based-gis-on-a-2-complex.md](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/development/ideas/2026-11-04-event-based-gis-on-a-2-complex.md)) reads the probability of a city's record of events as a contraction over its event complex:
- persistence on the cells (segments of worldlines, lenses of relations);
- kernels at the events;
- labels by modes of the states' dynamics, which are spins where the dynamics is symmetric.

**What the urban-science domain needed of physis.** Two general methods, neither naming a city:
- the contraction of such a record;
- a chain on a set classified at nested levels (divisions, groups, classes), whose modes say how long each level is remembered.

Neither existed. `master_equation.hpp` holds dense generators and their stationary laws. `lumping.hpp` holds the density of a fully observed path on one chain, not of many chains joined at events. `gauge/contraction.hpp` contracts Gibbs weights of ℤ₂ connections, not stochastic kernels along time.

## 2. Decision

### 2.1 Chains on a Classified Set (`statistical_mechanics/nested_chain.hpp`)

**The setting.**
- States are classified at levels ℓ = 0, …, L, each level refining the one before: the whole set at level 0, the states themselves at level L.
- π is a law on the states.
- At rate λ_ℓ a state is drawn again from π within its level-ℓ block.

**The theorem.** The generator L = Σ_ℓ λ_ℓ(Π_ℓ − I), with Π_ℓ the conditional expectation onto functions of the level-ℓ block, has these modes:
- V_0, the constants;
- V_m, the functions of the level-m block with mean zero on every level-(m−1) block.

L acts on V_m as −μ_m, with μ_m = λ_0 + … + λ_{m−1}. Hence:
- a level is never remembered longer than a coarser one;
- the propagator is closed: P_τ = Σ_m e^{−μ_m τ}(Π_m − Π_{m−1}).

This holds for any uneven tree. It is the design's Proposition 3.2 generalised, since a balanced tree with a uniform law is that proposition up to its parameters.

**What the header holds.**
- the classification, validated;
- the chain's memories, its kernel and its modes;
- the weight with which each mode carries a transition;
- the log-likelihood of observed transitions and its gradient;
- the rates' maximum likelihood, by Newton's method in log λ with Laplace's errors; a rate the data drive to zero is reported at its bound;
- the expected number of changes between two observations, by Hobolth and Jensen's endpoint conditioning (2005), whose integral the modes make closed.

### 2.2 The Probability of a Record (`statistical_mechanics/event_chains.hpp`)

**The model.**
- **Cells.** Each has a finite state space and a propagator over its life.
- **Pieces,** in the order of time. Each joins the cells ending there to those beginning there by a dense kernel: a law, a sum, a transition, an observation's likelihood, or a rule over several cells.

**The contraction.** It keeps the live cells in factors joined only where pieces join them, multiplies the factors of a piece's in-cells, sums the kernel against them, and carries each out-cell by its propagator. It rescales as it goes, refuses a factor beyond a limit, and returns log Z.

The design's Theorem 4.2 is this contraction, and its Theorem 4.5 (slicing time changes nothing) and Proposition 4.6 (Z = 1 without observations) are among the checks.

## 3. Alternatives Considered

- **A general tensor-network contraction with a found order** (decision 0007's plan). It is not needed while cells meet only at pieces. Along the order of time, the factors are the clusters events have joined, and the cost is their size.
- **The mode expansion inside the contraction.** The modes depend on each cell's generator. Where a model has them in closed form, as the nested chain does, the domain reads them there. A diagonalisation inside the contraction would serve none of the first uses.
- **A dense generator with calx's expm for the classified chain.** It is exact but costs N³ a time step, where the closed form costs L a pair, and the closed form shows why coarser levels outlast finer ones. calx's expm remains its check.

## 4. Consequences

- **Checks.** `statistical_mechanics` gains two headers and 16 checks:
  - `check_nested_chain.cpp`, 11: the closed form against expm, the modes as eigenspaces, the expected changes against Van Loan's integral, rates recovered from simulation, a level at its bound;
  - `check_event_chains.cpp`, 5: against enumeration, refinement, Z = 1, factors kept apart, faulty records refused.
- **Its first use.** The infrastructure's urban-science domain fits Melbourne's censuses of uses with the nested chain and contracts the buildings' chains on their event complex: the product of their transitions, exactly.
- **Planned.**
  - kernels that factor within a piece;
  - truncation of factors by weight;
  - forward filtering with backward sampling where factors grow too large;
  - the gauge models of relational states, with the history domain's order of time moved down beside them.

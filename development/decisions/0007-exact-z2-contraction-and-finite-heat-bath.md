# 0007. The Exact Contraction of ℤ₂ Gauge Theories, Couplings by 2-Cell, and the Heat Bath of Finite Groups

- **Status:** Accepted, provisionally (the contraction's order is the caller's, and its limits are set by hand)
- **Date:** 2026-10-24
- **Decided by:** the user, who asked that the spin foam of a history be computed exactly, its couplings fitted from the record, and then that the history domain reuse what the foundations already hold instead of writing its own: "Reusable, integrative science foundation" is the purpose of the repository. Claude moved the general methods here and wrote this record; the history domain's use of them is the infrastructure's [ADR 0042](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/development/decisions/0042-the-history-as-a-field-theory.md).

## 1. Context

The history domain's spin foam is a ℤ₂ lattice gauge theory on a 2-complex, the refined history complex K̂, with three features physis's `gauge` did not have:

- a coupling for each 2-cell (a tie's persistence by its colour, the balance of triangles), where `GaugeField` had one β;
- links fixed by observation (the attested signs) and by a gauge (the temporal edges), where every link was free;
- exact results on complexes far from surfaces, where physis knew Migdal and Witten's formula for closed surfaces and discs only.

The domain had written its own transfer matrix, its own Monte Carlo and its own SU(2) invariant dimensions. Each is general: none names a history. They belong in physis, where any 2-complex can use them.

## 2. Decision

### 2.1 Couplings by 2-Cell and Frozen Links

`GaugeField` takes a weight w_f for each 2-cell (`set_couplings`), so that the action is β Σ_f w_f (1 − Re tr U_f / N), and links may be frozen at a value (`freeze`, `thaw`, `frozen`). Metropolis skips frozen links. The heat bath of the staple groups applies where a link's 2-cells share one coupling (`link_coupling`), and otherwise falls back to Metropolis. Overrelaxation skips both kinds of link.

### 2.2 The Exact Contraction (`contraction.hpp`)

A `Z2Theory` is a set of links, the links each 2-cell's holonomy multiplies, a coupling per 2-cell and the fixed links. `Z2Theory::of(X, β)` builds it from a CW complex: a 2-cell multiplies the edges its word meets an odd number of times. `contract(T, order, …)` sums the free links exactly along an order of the 2-cells.

- **States** are kept in the character basis of ℤ₂^live (Drouffe and Zuber 1983, §2). A 2-cell multiplies by cosh β + c sinh β χ_V, and summing out a link keeps the characters without it.
- **Spin networks.** With no link fixed, the prior's states are confined to the ℤ/2 1-cycles of the cut graph, which are the boundary ℤ₂ spin networks. The cost is set by the widest cut, as for any transfer matrix (Baxter 1982, ch. 7).
- **Readings** come from forward and backward states paired at one point of the order: the mean of every link and 2-cell, two-point functions by an insertion carried along, and at chosen cuts the mean products of sets of live links. Where few links are live there, the law of how many products are negative is also given, along with the heaviest spin networks and their weights.
- **Derivatives.** `log_z_derivatives<M, K>` carries calx jets (`Taylor<M, K>`) through one forward pass. The gradient of log Z in M parameters is Σ ⟨h_f⟩ over each parameter's 2-cells, and the Hessian is their covariance, both exact.
- **Limits.** A state beyond `ContractionLimits` throws `Intractable`, and the caller samples instead.
- **Measure.** log Z is taken with the counting measure on the free links. With Haar measures of mass 1 (exact.hpp) it is less by E_free log 2.

### 2.3 The Heat Bath of Finite Groups

`heat_bath_sweep` for any `FiniteGaugeGroup` (ℤ_N, or a finite group by its table) draws U_e among the |G| elements with probability ∝ e^{β Σ_{f ∋ e} w_f Re tr U_f / N}. Each weight is computed by setting the element in turn, so the update is exact on any complex and for any couplings, including an edge met twice in a word.

It is also aperiodic. Metropolis on ℤ₂ has the flip as its only proposal, so a link on which the weight does not depend flips at every sweep, and a loop of an even number of such links never changes. `run()` uses this heat bath for finite groups.

### 2.4 Invariant Dimensions of SU(2)

`invariant_dimension(two_j)` in `symmetry/recoupling.hpp` gives dim Inv(V_{j₁} ⊗ … ⊗ V_{jₙ}), exactly as a morphe `Integer`, by coupling the spins in turn.

## 3. Alternatives Considered

- **Keeping the contraction in the history domain.** It names no object of history; a ℤ₂ theory on another 2-complex (a surface, a dual complex, a coloured graph) would have had to copy it. The domain keeps only what is the record's: which cells, which couplings, which links the record fixes, the order of time.
- **Fixing ℤ₂'s Metropolis instead of adding a heat bath.** A lazy proposal (stay with probability ½) would make Metropolis aperiodic, but at half the rate. The heat bath is exact for every finite group and needs no step size.
- **A general tensor-network contraction (any group, any order found automatically).** Wanted, and planned: other abelian groups in the character basis, and orders found by tree decomposition of the dual graph. ℤ₂ with the caller's order is what the first use needs, and it is checked exactly.

## 4. Consequences

- `gauge` gains `contraction.hpp`, couplings and frozen links in `field.hpp`, and the finite heat bath in `updates.hpp`. `group.hpp` gains the concept `FiniteGaugeGroup`, and `FiniteGroupAdapter` gains `order()`.
- `check_contraction.cpp` checks the contraction against brute force on a 3 × 3 torus with a coupling for each plaquette and three fixed links: log Z, every ⟨h⟩ and ⟨u⟩, 170 two-point functions and a cut's law. It also checks it against Migdal and Witten on the 4 × 4 torus and the surface of genus 2, and the flat connections at strong coupling. It checks the cycle confinement of the spin networks, the jets' gradient and Hessian, Metropolis and the heat bath with couplings and frozen links, and a circle in no 2-cell, which the heat bath finds and Metropolis cannot.
- `check_gauge.cpp` adds the heat bath for ℤ_3 and ℤ_2 against the character formula. `check_recoupling.cpp` adds invariant dimensions against the Catalan numbers and the triangle conditions.
- The history domain composes these; its checks compare its model in temporal gauge with this theory on the whole of K̂ on every canon within reach.

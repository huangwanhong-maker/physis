# quantum: Quantum Theory

*Status: draft (2026-10-23): the open two-level system.* Hilbert spaces, states and operators, spins and qubits in general, measurement, Liouvillians and mean-field closures are planned (the infrastructure's [map of the mathematics](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/docs/mathematics.md), §5).

## 1. Purpose

Quantum states and their dynamics, closed and open, for the sciences that need them: quantum optics, quantum information, chemistry, and the quantum neuromorphic models of neuroscience.

## 2. Concepts

- **An open system** is described by its density matrix, which obeys the GKSL master equation: a Hamiltonian part and, for each channel of loss, a jump operator L (Gorini, Kossakowski and Sudarshan 1976; Lindblad 1976). For a two-level system the density matrix is a Bloch vector in the unit ball, and the master equation is an affine flow on it, the optical Bloch equations, integrated exactly as a linear system of physis `dynamics`.
- **Quantum trajectories** unravel the same evolution into pure states (Dalibard, Castin and Mølmer 1992; Carmichael 1993). Between detected jumps the unnormalised state evolves by the non-Hermitian effective Hamiltonian, and its squared norm is the probability that no jump has yet occurred. At a jump the state is reset by L. Averaged over trajectories, it is the density matrix.
- **The two descriptions agree.** The mean waiting time between photons of a trajectory is the inverse of the master equation's emission rate, a bond the checks test directly.

## 3. Content

| Header | Content |
|---|---|
| `two_level.hpp` | `TwoLevel` (Omega, Delta, gamma, gamma_phi); `exp2`, the exponential of a 2×2 complex matrix by Cayley and Hamilton; `effective_hamiltonian`, `no_jump_propagator`; `Bloch`, the master equation as a linear system; `stationary`, `excited`, `emission_rate` |

## 4. Checks

`tests/check_quantum.cpp`: `exp2` against calx's `expm` of the real form, on random matrices and a Jordan block; the propagator between detections against its closed form on resonance; the stationary Bloch vector against the exact flow run long, and, without dephasing, against the textbook population; the Bloch ball kept; the mean waiting time against the inverse emission rate.

## 5. Dependencies

physis `dynamics`; calx `linalg` (through `dynamics`), and in the checks `quadrature` and `random`.

## 6. References

- Gorini, V., Kossakowski, A. and Sudarshan, E. C. G. (1976), Completely positive dynamical semigroups of N-level systems, *J. Math. Phys.* 17, 821–825.
- Lindblad, G. (1976), On the generators of quantum dynamical semigroups, *Commun. Math. Phys.* 48, 119–130.
- Dalibard, J., Castin, Y. and Mølmer, K. (1992), Wave-function approach to dissipative processes in quantum optics, *Phys. Rev. Lett.* 68, 580–583.
- Carmichael, H. (1993), *An Open Systems Approach to Quantum Optics*, Springer.
- Plenio, M. B. and Knight, P. L. (1998), The quantum-jump approach to dissipative dynamics in quantum optics, *Rev. Mod. Phys.* 70, 101–144.
- Cohen-Tannoudji, C. and Dalibard, J. (1986), Single-atom laser spectroscopy: looking for dark periods in fluorescence light, *Europhys. Lett.* 1, 441–448.
- Breuer, H.-P. and Petruccione, F. (2002), *The Theory of Open Quantum Systems*, Oxford University Press.

## 7. Open Questions

1. **N-level systems and many qubits:** states and operators of general dimension, tensor products, Liouvillians as matrices, and their sparse exponentials through calx.
2. **Spin models** (Ising, transverse-field) with open dynamics, shared with `statistical_mechanics`, and their mean-field closures.

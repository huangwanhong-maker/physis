// The defining constants of the SI, exact by definition since 20 May 2019 (BIPM 2019, The International System
// of Units, 9th ed., §2.2, Table 1). Measured constants (G, the fine-structure constant, masses) are CODATA's
// and versioned; they are not here (decision 0022).
#pragma once

namespace physis::units::si {

inline constexpr double caesium_frequency = 9192631770.0;  // Δν_Cs, Hz
inline constexpr double speed_of_light = 299792458.0;      // c, m/s
inline constexpr double planck = 6.62607015e-34;           // h, J s
inline constexpr double elementary_charge = 1.602176634e-19;  // e, C
inline constexpr double boltzmann = 1.380649e-23;          // k, J/K
inline constexpr double avogadro = 6.02214076e23;          // N_A, 1/mol
inline constexpr double luminous_efficacy = 683.0;         // K_cd, lm/W

}  // namespace physis::units::si

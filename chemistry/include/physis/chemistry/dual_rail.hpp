// Polynomial ODEs compiled to reaction networks of mass action, by dual rails (decision 0022).
//
// Any system x' = P(x), P polynomial, becomes a network of mass action whose concentrations are non-negative:
// each variable is the difference of two species, xᵢ = pᵢ − qᵢ (dual rails), and expanding every monomial of P in
// the pᵢ and qᵢ gives terms c · M(p, q) whose sign decides the rail they feed. A term c M with c > 0 in the
// equation of xᵢ is the reaction M → M + Pᵢ at rate c; with c < 0, M → M + Qᵢ at rate −c (the species of M are
// catalysts, unchanged); a constant term is ∅ → Pᵢ or ∅ → Qᵢ. Then d(pᵢ − qᵢ)/dt = Pᵢ(p − q) exactly, whatever
// the rails' common part does, and an annihilation Pᵢ + Qᵢ → ∅ at rate η keeps that part bounded without
// changing the difference (Fages, Le Guludec, Bournez and Pouly 2017, CMSB, LNCS 10545, 108; the encoding of
// signed values by pairs is older, e.g. Oishi and Klavins 2011, IET Syst. Biol. 5, 252).
//
// Terms with the same reactants and the same target variable are merged before their sign is read, so that the
// network has the fewest reactions this construction gives. Coefficients are rationals, so the network is exact.
// Checks: tests/check_dual_rail.cpp.
#pragma once

#include <cstddef>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <morphe/arithmetic/rational.hpp>

#include "kinetics.hpp"

namespace physis::chemistry {

struct Monomial {
  morphe::arithmetic::Rational coefficient = 0;
  std::vector<unsigned> exponents;  // one per variable
};

struct PolynomialSystem {
  std::vector<std::string> variables;
  std::vector<std::vector<Monomial>> rates;  // xᵢ' = Σ monomials of rates[i]

  double evaluate(std::size_t i, const std::vector<double>& x) const {
    double s = 0;
    for (const auto& m : rates[i]) {
      double v = m.coefficient.to_double();
      for (std::size_t j = 0; j < x.size(); j++)
        for (unsigned e = 0; e < m.exponents[j]; e++) v *= x[j];
      s += v;
    }
    return s;
  }
};

struct DualRailNetwork {
  Network network;
  std::vector<std::size_t> plus, minus;  // the species pᵢ and qᵢ of each variable
};

inline DualRailNetwork dual_rail(const PolynomialSystem& P, const morphe::arithmetic::Rational& annihilation = 1) {
  using morphe::arithmetic::Rational;
  const std::size_t n = P.variables.size();
  if (P.rates.size() != n) throw std::invalid_argument("physis::chemistry::dual_rail: one rate per variable");
  DualRailNetwork out;
  for (std::size_t i = 0; i < n; i++) out.plus.push_back(2 * i), out.minus.push_back(2 * i + 1);
  for (std::size_t i = 0; i < n; i++) out.network.species.push_back(P.variables[i] + "+"), out.network.species.push_back(P.variables[i] + "-");
  // (target variable, reactants as exponents of the 2n species) → net coefficient.
  std::map<std::pair<std::size_t, std::vector<unsigned>>, Rational> terms;
  for (std::size_t i = 0; i < n; i++)
    for (const auto& m : P.rates[i]) {
      if (m.exponents.size() != n) throw std::invalid_argument("physis::chemistry::dual_rail: one exponent per variable");
      // Expand Π_j (p_j − q_j)^{a_j}: choose k_j of the a_j factors to be −q_j.
      std::vector<unsigned> k(n, 0);
      while (true) {
        Rational c = m.coefficient;
        std::vector<unsigned> reactants(2 * n, 0);
        int sign = 1;
        for (std::size_t j = 0; j < n; j++) {
          // C(a_j, k_j)
          Rational binom = 1;
          for (unsigned t = 0; t < k[j]; t++) binom = binom * Rational(static_cast<long long>(m.exponents[j] - t)) / Rational(static_cast<long long>(t + 1));
          c = c * binom;
          reactants[2 * j] = m.exponents[j] - k[j], reactants[2 * j + 1] = k[j];
          if (k[j] % 2) sign = -sign;
        }
        if (sign < 0) c = -c;
        terms[{i, reactants}] = terms[{i, reactants}] + c;
        // The next choice of k.
        std::size_t j = 0;
        while (j < n && k[j] == m.exponents[j]) k[j] = 0, j++;
        if (j == n) break;
        k[j]++;
      }
    }
  for (const auto& [key, c] : terms) {
    const auto& [i, reactants] = key;
    if (c == Rational(0)) continue;
    Reaction r;
    const bool positive = Rational(0) < c;
    r.name = P.variables[i] + (positive ? "+ from " : "- from ");
    for (std::size_t s = 0; s < reactants.size(); s++)
      if (reactants[s]) r.reactants.push_back({s, reactants[s]}), r.products.push_back({s, reactants[s]}), r.name += out.network.species[s] + "^" + std::to_string(reactants[s]) + " ";
    const std::size_t target = positive ? out.plus[i] : out.minus[i];
    bool merged = false;
    for (auto& [s, b] : r.products)
      if (s == target) b += 1, merged = true;
    if (!merged) r.products.push_back({target, 1});
    r.rate = positive ? c : -c;
    out.network.reactions.push_back(std::move(r));
  }
  for (std::size_t i = 0; i < n; i++) {
    Reaction a;
    a.name = "annihilation of " + P.variables[i];
    a.reactants = {{out.plus[i], 1}, {out.minus[i], 1}};
    a.rate = annihilation;
    out.network.reactions.push_back(std::move(a));
  }
  return out;
}

// Lorenz's system with rational parameters, as a polynomial system.
inline PolynomialSystem lorenz_polynomials(const morphe::arithmetic::Rational& sigma, const morphe::arithmetic::Rational& r, const morphe::arithmetic::Rational& b) {
  PolynomialSystem P;
  P.variables = {"x", "y", "z"};
  P.rates = {
      {{sigma, {0, 1, 0}}, {-sigma, {1, 0, 0}}},                // σ(y − x)
      {{r, {1, 0, 0}}, {-1, {1, 0, 1}}, {-1, {0, 1, 0}}},       // r x − x z − y
      {{1, {1, 1, 0}}, {-b, {0, 0, 1}}},                        // x y − b z
  };
  return P;
}

}  // namespace physis::chemistry

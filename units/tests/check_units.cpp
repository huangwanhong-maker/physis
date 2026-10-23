// Checks of physis/units: quantities read with their units come out in SI with the right
// dimension; prefixes and whole symbols are told apart; derived units reduce to the base
// units as the SI brochure defines them; a quantity of the wrong dimension is refused.

#include <physis/units/quantity.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <stdexcept>
#include <string>

using namespace physis::units;
using calx::verification::relative;

int main() {
  calx::verification::Suite s;
  struct Case { const char* text; double value; Dimension dimension; };
  const Case cases[] = {
      {"10 ms", 0.01, dim::time},
      {"87.8 pA", 87.8e-12, dim::current},
      {"-65 mV", -0.065, dim::voltage},
      {"250 pF", 250e-12, dim::capacitance},
      {"120 mS/cm^2", 120e-3 / 1e-4, dim::conductance / dim::area},
      {"1 uF/cm^2", 1e-6 / 1e-4, dim::capacitance / dim::area},
      {"8 kHz", 8000, dim::frequency},
      {"1/s", 1, dim::frequency},
      {"2.5 MOhm", 2.5e6, dim::resistance},
      {"3 mM", 3, dim::concentration},
      {"1 mol/L", 1000, dim::concentration},
      {"9.81 m/s^2", 9.81, dim::acceleration},
      {"1 kg m^2 / s^3 / A", 1, dim::voltage},
      {"5 cm", 0.05, dim::length},
      {"2 cd", 2, dim::luminous_intensity},
      {"1 T", 1, dim::magnetic_field},
      {"3 min", 180, dim::time},
      {"0.5", 0.5, dim::one},
  };
  int good = 0;
  std::string bad;
  for (const auto& c : cases) {
    const Quantity q = parse_quantity(c.text);
    if (relative(q.value, c.value) < 1e-15 && q.dimension == c.dimension) good++;
    else bad += std::string(c.text) + " -> " + std::to_string(q.value) + " " + to_string(q.dimension) + "; ";
  }
  s.require("units: quantities with prefixes and compound units are read in SI with their dimensions",
            good == static_cast<int>(std::size(cases)), bad.empty() ? std::to_string(good) + " cases" : bad);
  s.require("units: the volt is kg m^2 s^-3 A^-1, the farad kg^-1 m^-2 s^4 A^2 (SI brochure, Table 4)",
            dim::voltage == Dimension(1, 2, -3, -1) && dim::capacitance == Dimension(-1, -2, 4, 2) && to_string(dim::conductance) == "kg^-1 m^-2 s^3 A^2");
  s.require("units: dimensions form a group (V = A Ohm, F S^-1 = s)", dim::voltage == dim::current * dim::resistance && dim::capacitance / dim::conductance == dim::time);
  bool refused = false;
  try { parse_as("10 nS", dim::current, "a current"); } catch (const std::invalid_argument&) { refused = true; }
  bool unknown = false;
  try { parse_quantity("3 furlongs"); } catch (const std::invalid_argument&) { unknown = true; }
  s.require("units: a conductance is refused where a current is meant, and an unknown unit is refused", refused && unknown);
  s.check("units: parse_as gives the SI value", parse_as("-70 mV", dim::voltage) + 0.070, 1e-17);
  return s.done();
}

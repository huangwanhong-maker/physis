// Physical dimensions and quantities.
//
// A dimension is an element of the free abelian group Z^7 over the SI base units (kg, m, s,
// A, K, mol, cd): units multiply by adding exponents (BIPM 2019, The International System
// of Units, 9th ed., §2.3.3). A quantity is a value in SI units with its dimension. Every
// value inside the libraries is in SI units; a quantity read from text ("10 ms", "87.8 pA",
// "120 mS/cm^2") is converted to SI at that boundary, and a parameter that declares its
// dimension refuses a quantity of another, so that a conductance cannot be given where a
// current is meant.
//
// Representation: dimensions are checked at run time, at the boundaries where quantities
// are read; compile-time dimensioned types may follow (physis's open decision 1).

#pragma once

#include <array>
#include <cctype>
#include <cmath>
#include <stdexcept>
#include <string>
#include <string_view>

namespace physis::units {

struct Dimension {
  // Exponents of kg, m, s, A, K, mol, cd.
  std::array<int, 7> e{};

  constexpr Dimension() = default;
  constexpr Dimension(int kg, int m, int s, int A, int K = 0, int mol = 0, int cd = 0) : e{kg, m, s, A, K, mol, cd} {}
  constexpr Dimension operator*(const Dimension& o) const {
    Dimension d;
    for (int i = 0; i < 7; i++) d.e[i] = e[i] + o.e[i];
    return d;
  }
  constexpr Dimension operator/(const Dimension& o) const {
    Dimension d;
    for (int i = 0; i < 7; i++) d.e[i] = e[i] - o.e[i];
    return d;
  }
  constexpr Dimension pow(int k) const {
    Dimension d;
    for (int i = 0; i < 7; i++) d.e[i] = e[i] * k;
    return d;
  }
  constexpr bool operator==(const Dimension&) const = default;
};

namespace dim {
inline constexpr Dimension one{};
inline constexpr Dimension mass{1, 0, 0, 0}, length{0, 1, 0, 0}, time{0, 0, 1, 0}, current{0, 0, 0, 1};
inline constexpr Dimension temperature{0, 0, 0, 0, 1}, amount{0, 0, 0, 0, 0, 1}, luminous_intensity{0, 0, 0, 0, 0, 0, 1};
inline constexpr Dimension area = length.pow(2), volume = length.pow(3);
inline constexpr Dimension frequency = one / time;
inline constexpr Dimension velocity = length / time, acceleration = velocity / time;
inline constexpr Dimension force = mass * acceleration, energy = force * length, power = energy / time;
inline constexpr Dimension pressure = force / area;
inline constexpr Dimension charge = current * time;
inline constexpr Dimension voltage = power / current;
inline constexpr Dimension capacitance = charge / voltage;
inline constexpr Dimension resistance = voltage / current, conductance = current / voltage;
inline constexpr Dimension magnetic_flux = voltage * time, magnetic_field = magnetic_flux / area, inductance = magnetic_flux / current;
inline constexpr Dimension concentration = amount / volume;
}  // namespace dim

// Multipliers for writing quantities in code: 10 * unit::ms.
namespace unit {
inline constexpr double s = 1, ms = 1e-3, us = 1e-6, min = 60, h = 3600;
inline constexpr double Hz = 1, kHz = 1e3;
inline constexpr double m = 1, cm = 1e-2, mm = 1e-3, um = 1e-6, cm2 = 1e-4, um2 = 1e-12;
inline constexpr double V = 1, mV = 1e-3, uV = 1e-6;
inline constexpr double A = 1, mA = 1e-3, uA = 1e-6, nA = 1e-9, pA = 1e-12;
inline constexpr double F = 1, uF = 1e-6, nF = 1e-9, pF = 1e-12;
inline constexpr double S = 1, mS = 1e-3, uS = 1e-6, nS = 1e-9;
inline constexpr double Ohm = 1, kOhm = 1e3, MOhm = 1e6, GOhm = 1e9;
inline constexpr double K = 1, mol = 1, mM = 1;  // mM: mmol/L = mol/m^3
}  // namespace unit

struct Quantity {
  double value = 0;  // SI
  Dimension dimension;
};

std::string to_string(const Dimension& d);
std::string symbol(const Dimension& d);
Quantity parse_quantity(std::string_view text);
double parse_as(std::string_view text, const Dimension& expected, std::string_view what = "");

// ---- implementation ----------------------------------------------------------------------------

namespace detail {

struct UnitEntry {
  std::string_view symbol;
  double scale;
  Dimension dimension;
  bool prefixable;
};

// Whole symbols are matched before prefixes, so that "m" is the metre, "min" the minute,
// "cd" the candela and "T" the tesla, and "mm", "TW", "kcd" are prefixed.
inline constexpr std::array<UnitEntry, 27> UNITS{{
    {"1", 1, dim::one, false},
    {"m", 1, dim::length, true},          {"g", 1e-3, dim::mass, true},       {"s", 1, dim::time, true},
    {"A", 1, dim::current, true},         {"K", 1, dim::temperature, true},   {"mol", 1, dim::amount, true},
    {"cd", 1, dim::luminous_intensity, true},
    {"Hz", 1, dim::frequency, true},      {"N", 1, dim::force, true},         {"Pa", 1, dim::pressure, true},
    {"J", 1, dim::energy, true},          {"W", 1, dim::power, true},         {"C", 1, dim::charge, true},
    {"V", 1, dim::voltage, true},         {"F", 1, dim::capacitance, true},   {"Ohm", 1, dim::resistance, true},
    {"ohm", 1, dim::resistance, true},    {"\xCE\xA9", 1, dim::resistance, true},
    {"S", 1, dim::conductance, true},     {"Wb", 1, dim::magnetic_flux, true}, {"T", 1, dim::magnetic_field, true},
    {"H", 1, dim::inductance, true},      {"L", 1e-3, dim::volume, true},     {"M", 1e3, dim::concentration, true},
    {"min", 60, dim::time, false},        {"h", 3600, dim::time, false},
}};

inline constexpr std::array<std::pair<std::string_view, double>, 13> PREFIXES{{
    {"T", 1e12}, {"G", 1e9}, {"M", 1e6}, {"k", 1e3}, {"d", 1e-1}, {"c", 1e-2}, {"m", 1e-3},
    {"u", 1e-6}, {"\xC2\xB5", 1e-6}, {"\xCE\xBC", 1e-6}, {"n", 1e-9}, {"p", 1e-12}, {"f", 1e-15},
}};

inline bool lookup(std::string_view sym, double& scale, Dimension& d) {
  for (const auto& u : UNITS)
    if (u.symbol == sym) { scale = u.scale; d = u.dimension; return true; }
  for (const auto& [p, f] : PREFIXES) {
    if (sym.size() <= p.size() || sym.substr(0, p.size()) != p) continue;
    for (const auto& u : UNITS)
      if (u.prefixable && u.symbol == sym.substr(p.size())) { scale = f * u.scale; d = u.dimension; return true; }
  }
  return false;
}

}  // namespace detail

inline std::string to_string(const Dimension& d) {
  static constexpr const char* names[7] = {"kg", "m", "s", "A", "K", "mol", "cd"};
  std::string out;
  for (int i = 0; i < 7; i++) {
    if (!d.e[i]) continue;
    if (!out.empty()) out += " ";
    out += names[i];
    if (d.e[i] != 1) out += "^" + std::to_string(d.e[i]);
  }
  return out.empty() ? "1" : out;
}

// "<number> [unit]": a unit is a product of prefixed units with integer powers, and after a
// '/' all further units divide: "mS/cm^2", "V/s", "1/s", "mol/L", "kg m^2 / s^3".
inline Quantity parse_quantity(std::string_view text) {
  auto fail = [&](const char* why) { throw std::invalid_argument(std::string("quantity \"") + std::string(text) + "\": " + why); };
  std::size_t i = 0;
  while (i < text.size() && std::isspace(static_cast<unsigned char>(text[i]))) i++;
  const std::string str(text.substr(i));
  std::size_t used = 0;
  double value = 0;
  try { value = std::stod(str, &used); } catch (...) { fail("no number"); }
  if (!std::isfinite(value)) fail("not finite");
  std::string_view rest = std::string_view(str).substr(used);
  auto trim = [&] {
    while (!rest.empty() && std::isspace(static_cast<unsigned char>(rest.front()))) rest.remove_prefix(1);
    while (!rest.empty() && std::isspace(static_cast<unsigned char>(rest.back()))) rest.remove_suffix(1);
  };
  trim();
  Quantity q{value, dim::one};
  int sign = 1;
  if (!rest.empty() && rest.front() == '/') {
    sign = -1;
    rest.remove_prefix(1);
    trim();
    if (rest.empty()) fail("nothing after '/'");
  }
  while (!rest.empty()) {
    std::size_t j = 0;
    while (j < rest.size() && rest[j] != '/' && rest[j] != '*' && rest[j] != ' ' && rest[j] != '^' && rest[j] != '.') j++;
    const std::string_view sym = rest.substr(0, j);
    rest.remove_prefix(j);
    int power = 1;
    if (!rest.empty() && rest.front() == '^') {
      rest.remove_prefix(1);
      std::size_t k = 0;
      if (k < rest.size() && (rest[k] == '-' || rest[k] == '+')) k++;
      const std::size_t digits = k;
      while (k < rest.size() && std::isdigit(static_cast<unsigned char>(rest[k]))) k++;
      if (k == digits) fail("a power without digits");
      power = std::stoi(std::string(rest.substr(0, k)));
      rest.remove_prefix(k);
    }
    double scale = 1;
    Dimension d;
    if (sym.empty() || !detail::lookup(sym, scale, d)) fail("unknown unit");
    q.value *= std::pow(scale, sign * power);
    q.dimension = q.dimension * d.pow(sign * power);
    while (!rest.empty() && rest.front() == ' ') rest.remove_prefix(1);
    if (rest.empty()) break;
    if (rest.front() == '/') {
      sign = -1;
      rest.remove_prefix(1);
    } else if (rest.front() == '*' || rest.front() == '.') {
      rest.remove_prefix(1);
    }  // else: a space separated product, "kg m^2"
    while (!rest.empty() && rest.front() == ' ') rest.remove_prefix(1);
    if (rest.empty()) fail("a unit is missing at the end");
  }
  return q;
}

// The conventional symbol of a dimension's coherent SI unit where it has one (V, A, F, S,
// Hz, V/s, ...), else its product of base units.
inline std::string symbol(const Dimension& d) {
  struct Named { Dimension d; const char* s; };
  static const Named named[] = {
      {dim::one, "1"},         {dim::time, "s"},           {dim::frequency, "Hz"},      {dim::length, "m"},
      {dim::mass, "kg"},       {dim::current, "A"},        {dim::voltage, "V"},         {dim::capacitance, "F"},
      {dim::conductance, "S"}, {dim::resistance, "Ohm"},   {dim::charge, "C"},          {dim::energy, "J"},
      {dim::power, "W"},       {dim::pressure, "Pa"},      {dim::force, "N"},           {dim::temperature, "K"},
      {dim::amount, "mol"},    {dim::concentration, "mol/m^3"}, {dim::magnetic_field, "T"}, {dim::area, "m^2"},
      {dim::volume, "m^3"},    {dim::voltage / dim::time, "V/s"}, {dim::current / dim::time, "A/s"},
      {dim::conductance / dim::area, "S/m^2"}, {dim::capacitance / dim::area, "F/m^2"}, {dim::current / dim::area, "A/m^2"},
  };
  for (const auto& n : named) if (n.d == d) return n.s;
  return to_string(d);
}

inline double parse_as(std::string_view text, const Dimension& expected, std::string_view what) {
  const Quantity q = parse_quantity(text);
  if (!(q.dimension == expected))
    throw std::invalid_argument(std::string(what.empty() ? "a quantity" : what) + ": \"" + std::string(text) + "\" is in " +
                                symbol(q.dimension) + ", not " + symbol(expected) + " (dimension " + to_string(q.dimension) + ", not " + to_string(expected) + ")");
  return q.value;
}

}  // namespace physis::units

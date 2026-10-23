# units: Dimensions and Quantities

*Status: draft (2026-10-23): dimensions and quantities.* CODATA constants, versioned, and compile-time dimensioned types are planned.

## 1. Purpose

Every quantity at a public interface carries its units, and inside the libraries every value is in SI units (invariant 3). This module gives dimensions their algebra and reads quantities written with their units, so that the boundary where a number enters is also where its dimension is checked.

## 2. Concepts

The dimensions form the free abelian group ℤ⁷ over the SI base units kg, m, s, A, K, mol and cd; a unit is a positive multiple of a dimension, and units multiply by adding exponents (BIPM 2019, §2.3.3). A quantity is a value in SI units with its dimension.

## 3. Content

`include/physis/units/quantity.hpp`:

- `Dimension`, with products, quotients and integer powers; the named dimensions of `dim::` (voltage, capacitance, conductance, concentration, ...); `to_string`.
- `unit::` multipliers for writing quantities in code (`10 * unit::ms`).
- `parse_quantity("120 mS/cm^2")`: a number and a unit written as products and quotients of prefixed SI units with integer powers; whole symbols are matched before prefixes, so that `m` is the metre, `min` the minute, `cd` the candela, `T` the tesla and `M` the molar; `parse_as(text, dimension)` refuses a quantity of another dimension.

## 4. Checks

`tests/check_units.cpp`: eighteen quantities read into SI with their dimensions; the derived units against the SI brochure's Table 4; the group law; the refusal of wrong dimensions and unknown units.

## 5. Dependencies

None.

## 6. References

- Bureau International des Poids et Mesures (2019), *The International System of Units (SI)*, 9th ed.

## 7. Open Questions

1. Compile-time dimensioned types beside the run-time checks, and calx's typed IR annotations (physis's open decision 1).
2. Non-multiplicative units (degrees Celsius, decibels) at the boundaries.

#!/usr/bin/env python3

from sympy import dsolve, solve, symbols, Eq
from symplyphysics import units, convert_to, Quantity, print_expression
from symplyphysics.laws.hydro import inner_pressure_is_sum_of_pressures as inner_pressure
from symplyphysics.laws.hydro import hydrostatic_pressure_via_density_and_height as hydrostatic_pressure
from symplyphysics.laws.hydro import inner_pressure_of_fluid_is_constant as constant_pressure_law

# Description
## A U-tube contains two liquids at static equilibrium:
## 1) water of density rho_water = 998 kg/m^3 in the right arm and
## 2) oil of unknown density rho_oil in the left.
## The distance between the oil-water interface and the end of the water column is 135 mm.
## The oil column is 12.3 mm higher than the water column.
## Find the density of the oil.

atmospheric_pressure = symbols("p_atm")
water_density, oil_density = symbols("rho_water rho_oil")
water_column_height = symbols("h_water")
oil_column_difference = symbols("d_oil")

oil_column_height = water_column_height + oil_column_difference

given_values = {
    water_density: Quantity(998 * units.kilogram / units.meter**3),
    water_column_height: Quantity(135 * units.millimeter),
    oil_column_difference: Quantity(12.3 * units.millimeter),
}

# Left arm

hydrostatic_pressure_left_arm_law = hydrostatic_pressure.law.subs({
    hydrostatic_pressure.density: oil_density,
    hydrostatic_pressure.height: oil_column_height,
})
hydrostatic_pressure_left_arm = solve(hydrostatic_pressure_left_arm_law,
    hydrostatic_pressure.hydrostatic_pressure)[0]

inner_pressure_left_arm_law = inner_pressure.law.subs({
    inner_pressure.static_pressure: atmospheric_pressure,
    inner_pressure.dynamic_pressure: 0,
    inner_pressure.hydrostatic_pressure: hydrostatic_pressure_left_arm,
})
inner_pressure_left_arm = solve(inner_pressure_left_arm_law, inner_pressure.inner_pressure)[0]

# Right arm

hydrostatic_pressure_right_arm_law = hydrostatic_pressure.law.subs({
    hydrostatic_pressure.density: water_density,
    hydrostatic_pressure.height: water_column_height,
})
hydrostatic_pressure_right_arm = solve(hydrostatic_pressure_right_arm_law,
    hydrostatic_pressure.hydrostatic_pressure)[0]

inner_pressure_right_arm_law = inner_pressure.law.subs({
    inner_pressure.static_pressure: atmospheric_pressure,
    inner_pressure.dynamic_pressure: 0,
    inner_pressure.hydrostatic_pressure: hydrostatic_pressure_right_arm,
})
inner_pressure_right_arm = solve(inner_pressure_right_arm_law, inner_pressure.inner_pressure)[0]

# Use Bernoulli's principle to equate the total pressure in both arms
dsolved = dsolve(constant_pressure_law.law,
    constant_pressure_law.inner_pressure(constant_pressure_law.time))
dsolved_left = dsolved.subs(constant_pressure_law.inner_pressure(constant_pressure_law.time),
    inner_pressure_left_arm)
dsolved_right = dsolved.subs(constant_pressure_law.inner_pressure(constant_pressure_law.time),
    inner_pressure_right_arm)
solved_left = solve([dsolved_left, dsolved_right], (inner_pressure_left_arm, "C1"),
    dict=True)[0][inner_pressure_left_arm]
constant_pressure_eq = Eq(inner_pressure_left_arm, solved_left)

oil_density_formula = solve(constant_pressure_eq, oil_density)[0]

print(f"Formula for oil density:\n{print_expression(Eq(oil_density, oil_density_formula))}")

# Note that the formula for oil density does not depend on either atmospheric pressure or
# the free-fall acceleration

assert atmospheric_pressure not in oil_density_formula.free_symbols
assert units.acceleration_due_to_gravity not in oil_density_formula.free_symbols

oil_density_quantity = Quantity(oil_density_formula.subs(given_values))
oil_density_value = convert_to(oil_density_quantity, units.kilogram / units.meter**3).evalf(3)

print(f"Oil density is {oil_density_value} kg/m^3")

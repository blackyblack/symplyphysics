#!/usr/bin/env python3

from sympy import solve, symbols, dsolve, Eq
from symplyphysics import units, convert_to, prefixes, Quantity, print_expression
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.hydro import inner_pressure_of_fluid as inner_pressure
from symplyphysics.laws.hydro import inner_pressure_of_fluid_is_constant as bernoullis_equation
from symplyphysics.laws.hydro import hydrostatic_pressure_from_density_and_depth as hydrostatic_pressure

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
    hydrostatic_pressure.depth: oil_column_height,
})
hydrostatic_pressure_left_arm = solve(
    hydrostatic_pressure_left_arm_law,
    hydrostatic_pressure.hydrostatic_pressure
)[0]

inner_pressure_left_arm_law = inner_pressure.law.subs({
    inner_pressure.static_pressure: atmospheric_pressure,
    inner_pressure.dynamic_pressure: 0,
    inner_pressure.hydrostatic_pressure: hydrostatic_pressure_left_arm,
})
inner_pressure_left_arm = solve(inner_pressure_left_arm_law, inner_pressure.inner_pressure)[0]

# Right arm

hydrostatic_pressure_right_arm_law = hydrostatic_pressure.law.subs({
    hydrostatic_pressure.density: water_density,
    hydrostatic_pressure.depth: water_column_height,
})
hydrostatic_pressure_right_arm = solve(
    hydrostatic_pressure_right_arm_law,
    hydrostatic_pressure.hydrostatic_pressure
)[0]

inner_pressure_right_arm_law = inner_pressure.law.subs({
    inner_pressure.static_pressure: atmospheric_pressure,
    inner_pressure.dynamic_pressure: 0,
    inner_pressure.hydrostatic_pressure: hydrostatic_pressure_right_arm,
})
inner_pressure_right_arm = solve(inner_pressure_right_arm_law, inner_pressure.inner_pressure)[0]

# The inner pressure in conserved at all points of the fluid

inner_pressure_solved = dsolve(
    bernoullis_equation.law, 
    bernoullis_equation.inner_pressure(bernoullis_equation.time)
)
assert expr_equals(inner_pressure_solved.rhs, symbols("C1"))

# Since the total pressure in conserved, we can equate the total pressure in both arms

equation = Eq(inner_pressure_left_arm, inner_pressure_right_arm)
oil_density_formula = solve(equation, oil_density)[0]

print(f"Formula for oil density:\n{print_expression(Eq(oil_density, oil_density_formula))}")

# Note that the formula for oil density does not depend on either atmospheric pressure or
# the free-fall acceleration

assert atmospheric_pressure not in oil_density_formula.free_symbols
assert units.acceleration_due_to_gravity not in oil_density_formula.free_symbols

oil_density_quantity = Quantity(oil_density_formula.subs(given_values))
oil_density_value = convert_to(oil_density_quantity, units.kilogram / units.meter**3).evalf(3)

print(f"Oil density is {oil_density_value} kg/m^3")

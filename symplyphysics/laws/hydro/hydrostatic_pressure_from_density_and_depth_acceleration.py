from sympy import Eq, solve
from symplyphysics import units, Quantity, Symbol, print_expression, validate_input, validate_output

# Description
# Law: P = ρ * a * h
# P - hydrostatic pressure
# ρ - density of liquid
# a - acceleration of vessel (should be directed vertically upwards)
# h - depth

# Conditions
# - No other accelerations involved, ie gravitational acceleration

density = Symbol("density", units.mass / units.volume)
depth = Symbol("depth", units.length)
acceleration = Symbol("acceleration", units.acceleration)
hydrostatic_pressure = Symbol("hydrostatic_pressure", units.pressure)

law = Eq(hydrostatic_pressure, density * acceleration * depth)


def print_law() -> str:
    return print_expression(law)


@validate_input(density_=density, depth_=depth, acceleration_=acceleration)
@validate_output(hydrostatic_pressure)
def calculate_hydrostatic_pressure(density_: Quantity, depth_: Quantity, acceleration_: Quantity) -> Quantity:
    result_pressure_expr = solve(law, hydrostatic_pressure, dict=True)[0][hydrostatic_pressure]
    result_expr = result_pressure_expr.subs({density: density_, depth: depth_, acceleration: acceleration_,})
    return Quantity(result_expr)

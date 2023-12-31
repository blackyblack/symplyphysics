from sympy import Eq, solve
from symplyphysics import units, Quantity, Symbol, print_expression, validate_input, validate_output

# Description
# Law: P = ρ * g * h
# P - hydrostatic pressure
# ρ - density of liquid
# g - acceleration of gravity
# h - depth

density = Symbol("density", units.mass / units.volume)
depth = Symbol("depth", units.length)
hydrostatic_pressure = Symbol("hydrostatic_pressure", units.pressure)

law = Eq(hydrostatic_pressure, density * units.acceleration_due_to_gravity * depth)


def print_law():
    return print_expression(law)


@validate_input(density_=density, depth_=depth)
@validate_output(hydrostatic_pressure)
def calculate_hydrostatic_pressure(density_: Quantity, depth_: Quantity) -> Quantity:
    result_pressure_expr = solve(law, hydrostatic_pressure, dict=True)[0][hydrostatic_pressure]
    result_expr = result_pressure_expr.subs({density: density_, depth: depth_})
    return Quantity(result_expr)

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)

# Description
## When liquid flows, it causes additional pressure, known as dynamic pressure.
## Law: p = ro * V**2 / 2, where
## ro is density of liquid,
## V is flow velocity.

# Conditions
## Liquid is ideal: no heat losses, no liquid friction losses and liquid is not compressible.

liquid_density = Symbol("liquid_density", units.mass / units.volume)
flow_velocity = Symbol("flow_velocity", units.velocity)
dynamic_pressure = Symbol("dynamic_pressure", units.pressure)

law = Eq(dynamic_pressure, liquid_density * flow_velocity**2 / 2)


def print_law() -> str:
    return print_expression(law)


@validate_input(density_=liquid_density, velocity_=flow_velocity)
@validate_output(dynamic_pressure)
def calculate_pressure(density_: Quantity, velocity_: Quantity) -> Quantity:
    result_pressure_expr = solve(law, dynamic_pressure, dict=True)[0][dynamic_pressure]
    result_expr = result_pressure_expr.subs({liquid_density: density_, flow_velocity: velocity_})
    return Quantity(result_expr)

from math import pi

from sympy import Eq, solve

from symplyphysics import (Quantity, Symbol, print_expression, units,
                           validate_input, validate_output)

# Law: delta_p = 8 * mu * L * Q / (pi * R**4)
# delta_p - pressure difference,
# mu is viscosity,
# L is length,
# Q is the volumetric flow rate,
# R is the pipe radius
# A is cross-section area of pipe.

dynamic_viscosity = Symbol("dynamic_viscosity", units.pressure * units.time)
length = Symbol("length", units.length)
flow_rate = Symbol("flow_rate", units.volume / units.time)
radius = Symbol("radius", units.length)
delta_pressure = Symbol("delta_pressure", units.pressure)


law = Eq(delta_pressure, 8 * dynamic_viscosity *
         length * flow_rate / (pi * radius**4))


def print_law() -> str:
    return print_expression(law)


@validate_input(
    dynamic_viscosity_=dynamic_viscosity,
    length_=length,
    flow_rate_=flow_rate,
    radius_=radius,
)
@validate_output(delta_pressure)
def calculate_delta_p(
        dynamic_viscosity_: Quantity,
        length_: Quantity,
        flow_rate_: Quantity,
        radius_: Quantity) -> Quantity:
    result_expr = solve(law, delta_pressure, dict=True)[0][delta_pressure]
    result_applied = result_expr.subs({
        dynamic_viscosity: dynamic_viscosity_,
        length: length_,
        flow_rate: flow_rate_,
        radius: radius_,
    })
    return Quantity(result_applied)

from sympy import (Eq, solve)
from symplyphysics import (
    units, Quantity, Symbol, print_expression, validate_input, validate_output, dimensionless, convert_to)


# Description
# The Reynolds number is a dimensionless quantity that characterizes the flow of a fluid in a pipe.
# Law: Re = rho * d * v / mu, where
# rho is fluid density
# d is hydraulic diameter of the pipe
# v is velocity of the fluid.


diameter = Symbol("diameter", units.length)
density = Symbol(
    "density", (units.mass / units.volume)
)
velocity = Symbol("velocity", units.velocity)
dynamic_viscosity = Symbol(
    "dynamic_viscosity", units.pressure * units.time
)

reynolds_number = Symbol("reynolds_number", dimensionless)

law = Eq(reynolds_number, density * velocity * diameter / dynamic_viscosity)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    diameter_=diameter,
    density_=density,
    velocity_=velocity,
    dynamic_viscosity_=dynamic_viscosity
)
@validate_output(reynolds_number)
def calculate_reynolds_number(diameter_: Quantity, density_: Quantity,
                              velocity_: Quantity, dynamic_viscosity_: Quantity) -> float:
    result_expr = solve(law, reynolds_number, dict=True)[0][reynolds_number]
    result_applied = result_expr.subs({
        diameter: diameter_,
        density: density_,
        velocity: velocity_,
        dynamic_viscosity: dynamic_viscosity_
    })
    result = Quantity(result_applied)
    
    return float(convert_to(result, Quantity(dimensionless)))

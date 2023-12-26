from sympy import (Eq, solve)
from symplyphysics import (
    units, Quantity, Symbol, print_expression, validate_input, validate_output, dimensionless)


# Description
# This file contains the implementation of the Reynolds number calculation in fluid dynamics.
## The Reynolds number quantifies the relative importance of these two types of forces for given flow conditions,
## and is a guide to when turbulent flow will occur in a particular situation
# The Reynolds number is a dimensionless quantity that characterizes the flow of a fluid in a pipe.
# It is calculated using the hydraulic diameter, fluid density, fluid velocity, and dynamic viscosity.
# Law: Re = rho * d * v / mu, where
# rho is fluid density
# d is hydraulic diameter of the pipe
# v is velocity of the fluid m/s.


diameter = Symbol("diameter", units.length)  # hydraulic diameter of the pipe
density = Symbol(
    "density", (units.mass / units.volume)  # density of the fluid kg/m^3
)
velocity = Symbol("volume", units.velocity)  # velocity of the fluid m/s
dynamic_viscosity = Symbol(
    # dynamic viscosity of the fluid Pa*s
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
    return Quantity(result_applied)

"""

"""

from sympy import Eq, sqrt
from sympy.physics.units import speed_of_light
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## The Lorentz factor, a.k.a. Lorentz term or gamma factor, is a quantity that expresses
## how much the measurements of time, length, and other physical properties change for
## an object while that object is moving.

# Definition: gamma = 1 / sqrt(1 - v**2/c**2)
## gamma - Lorentz factor
## v - object's velocity
## c - speed of light

lorentz_factor = Symbol("lorentz_factor", dimensionless)
velocity = Symbol("velocity", units.velocity)

definition = Eq(lorentz_factor, 1 / sqrt(1 - velocity**2 / speed_of_light**2))


def print_law() -> str:
    return print_expression(definition)


@validate_input(velocity_=velocity)
@validate_output(lorentz_factor)
def calculate_lorentz_factor(velocity_: Quantity) -> Quantity:
    result = definition.rhs.subs(velocity, velocity_)
    return Quantity(result)

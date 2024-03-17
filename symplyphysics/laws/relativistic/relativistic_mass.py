from sympy import Eq, solve, sqrt
from sympy.physics.units import speed_of_light

from symplyphysics import (Quantity, Symbol, print_expression, units, validate_input,
    validate_output)

# Description
## The relativistic mass is the sum total quantity of energy in a body or system
## Law: m_rel = m / sqrt(1 - v**2 / c**2), where
## m_rel is relativistic mass,
## m is rest mass,
## v is velocity,
## c is speed of light.

# Conditions
## Non-zero rest mass

# This is equivalent of symbols.basic.mass
rest_mass = Symbol("rest_mass", units.mass)
velocity = Symbol("velocity", units.velocity)
relativistic_mass = Symbol("relativistic_mass", units.mass)

law = Eq(relativistic_mass, rest_mass / sqrt(1 - velocity**2 / speed_of_light**2))


def print_law() -> str:
    return print_expression(law)


@validate_input(rest_mass_=rest_mass, velocity_=velocity)
@validate_output(relativistic_mass)
def calculate_relativistic_mass(rest_mass_: Quantity, velocity_: Quantity) -> Quantity:
    result_expr = solve(law, relativistic_mass, dict=True)[0][relativistic_mass]
    mass_applied = result_expr.subs({rest_mass: rest_mass_, velocity: velocity_})
    return Quantity(mass_applied)

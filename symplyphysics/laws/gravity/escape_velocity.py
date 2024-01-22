from sympy import Eq, solve, sqrt
from sympy.physics.units import gravitational_constant
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)


# Law: V = √(G * M / (R + h))
# Where:
# V - initial velocity
# M - planet_mass of planet
# G - gravitational constant
# h - height above the planet surface
# R - radius of planet


velocity = Symbol("initial_velocity", units.velocity)
planet_mass = Symbol("planet_mass", units.mass)
radius = Symbol("radius", units.length)
height = Symbol("height", units.length)

law = Eq(velocity, sqrt(gravitational_constant * planet_mass / (radius + height)))


def print_law():
    print_expression(law)


@validate_input(
        planet_mass_=planet_mass,
        radius_=radius,
        height_=height)
@validate_output(velocity)
def calculate_velocity(planet_mass_: Quantity, radius_: Quantity, height_: Quantity) -> Quantity:
    result_velocity_expr = solve(law, velocity, dict=True)[0][velocity]
    result_expr = result_velocity_expr.subs({planet_mass: planet_mass_, radius: radius_, height: height_})
    return Quantity(result_expr)

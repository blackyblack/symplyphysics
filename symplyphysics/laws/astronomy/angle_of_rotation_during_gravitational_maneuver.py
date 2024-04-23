from sympy import Eq, solve, atan
from sympy.physics.units import gravitational_constant
from symplyphysics import (clone_symbol, symbols, units, Quantity, Symbol, print_expression,
    validate_input, validate_output, angle_type)

## Description
## A gravitational maneuver is a purposeful change in the trajectory and flight speed of a spacecraft under the influence
## of the gravitational fields of celestial bodies.
## The angle of the gravitational maneuver depends on the aiming range, the mass of the planet and the velocity of the rocket relative to the planet.
## The aiming range is the distance between the asymptote of the hyperbolic trajectory of the circumnavigation of the planet and its focus
## coinciding with the center of the planet.
## https://ru.wikipedia.org/wiki/Гравитационный_манёвр#/media/Файл:Gravity_assist_-_ru.svg

## Law is: phi = 2 * arctg(G * M / (c * v^2)), where
## phi - angle of rotation during a gravitational maneuver (angle at which the velocity vector of the rocket rotates),
## G - gravitational constant,
## M - mass of the planet,
## c - aiming range,
## v - rocket's velocity relative to the planet.

angle = Symbol("angle", angle_type)
planet_mass = clone_symbol(symbols.basic.mass, "planet_mass")
aiming_range = Symbol("aiming_range", units.length)
rocket_speed = Symbol("rocket_speed", units.velocity)

law = Eq(angle, 2 * atan(gravitational_constant * planet_mass / (aiming_range * rocket_speed**2)))


def print_law() -> str:
    return print_expression(law)


@validate_input(planet_mass_=planet_mass, aiming_range_=aiming_range, rocket_speed_=rocket_speed)
@validate_output(angle)
def calculate_angle(planet_mass_: Quantity, aiming_range_: Quantity,
    rocket_speed_: Quantity) -> Quantity:
    result_velocity_expr = solve(law, angle, dict=True)[0][angle]
    result_expr = result_velocity_expr.subs({
        planet_mass: planet_mass_,
        aiming_range: aiming_range_,
        rocket_speed: rocket_speed_
    })
    return Quantity(result_expr)

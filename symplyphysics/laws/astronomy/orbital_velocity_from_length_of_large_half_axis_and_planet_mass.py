from sympy import Eq, solve, sqrt
from sympy.physics.units import gravitational_constant
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, symbols, clone_symbol)

## Description
## The orbital velocity of a body (usually a planet, a natural or artificial satellite, or a multiple star) is the speed at which it rotates around the barycenter of the system, usually around a more massive body.
## Law: V = √(G * M * ((2 / r) - 1 / a))
## Where:
## V is orbital velocity (velocity of rotating body)
## M is mass of planet (mass of central body)
## G is gravitational constant
## r is the distance between the rotating body and the central body
## a - is length of the large half-axis

## Conditions
## r <= a

orbital_velocity = Symbol("orbital_velocity", units.velocity)
distance = Symbol("distance", units.length)
large_half_axis_length = Symbol("large_half_axis_length", units.length)
planet_mass = clone_symbol(symbols.basic.mass)

law = Eq(
    orbital_velocity,
    sqrt(gravitational_constant * planet_mass * ((2 / distance) - (1 / large_half_axis_length))))


def print_law() -> str:
    return print_expression(law)


@validate_input(planet_mass_=planet_mass,
    distance_=distance,
    large_half_axis_length_=large_half_axis_length)
@validate_output(orbital_velocity)
def calculate_orbital_velocity(planet_mass_: Quantity, distance_: Quantity,
    large_half_axis_length_: Quantity) -> Quantity:
    if distance_.scale_factor > large_half_axis_length_.scale_factor:
        raise ValueError(
            "The distance between the rotating body and the central body must be less or equal to the large half-axis."
        )
    result_velocity_expr = solve(law, orbital_velocity, dict=True)[0][orbital_velocity]
    result_expr = result_velocity_expr.subs({
        planet_mass: planet_mass_,
        distance: distance_,
        large_half_axis_length: large_half_axis_length_
    })
    return Quantity(result_expr)

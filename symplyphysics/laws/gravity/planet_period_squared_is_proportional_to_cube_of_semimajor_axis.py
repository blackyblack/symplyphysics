from sympy import Eq, solve, pi
from sympy.physics.units import gravitational_constant
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    clone_symbol,
    symbols,
)

# Description
## Also known as Kepler's third law of planetary motion, the law of periods relates
## the period of rotation of any planet to the semimajor axis of its orbit.

# Law: T**2 = 4*pi**2 / (G*M) * a**3
## T - planet's period of rotation
## G - gravitational constant
## M - mass of the attracting body, the Sun in case of the solar system
## a - semimajor axis of the planet's orbit, its radius in case of a round orbit

rotation_period = Symbol("rotation_period", units.time, positive=True)
attracting_mass = clone_symbol(symbols.basic.mass, "attracting_mass", positive=True)
semimajor_axis = Symbol("semimajor_axis", units.length, positive=True)

law = Eq(
    rotation_period**2,
    4 * pi**2 / (gravitational_constant * attracting_mass) * semimajor_axis**3,
)

# TODO: derive law from Newton's second law of motion


def print_law() -> str:
    return print_expression(law)


@validate_input(attracting_mass_=attracting_mass, semimajor_axis_=semimajor_axis)
@validate_output(rotation_period)
def calculate_rotation_period(
    attracting_mass_: Quantity,
    semimajor_axis_: Quantity,
) -> Quantity:
    result = solve(law, rotation_period)[0].subs({
        attracting_mass: attracting_mass_,
        semimajor_axis: semimajor_axis_,
    })
    return Quantity(result)

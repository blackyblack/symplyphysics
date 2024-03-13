from sympy import (Eq, solve,)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    dimensionless,
)

# Description
## The Titius—Bode rule is an empirical formula that roughly describes the distances between the planets of the
## Solar System and the Sun (the average radii of the orbits).

## Law is: R = 0.4 + 0.3 * 2^i, where
## R - radius of the planet of the solar system in astronomical units,
## i - number of planet.

radius_of_orbit = Symbol("radius_of_orbit", units.length)

number_of_planet = Symbol("number_of_planet", dimensionless)

first_constant = Quantity(0.4 * units.astronomical_unit)
second_constant = Quantity(0.3 * units.astronomical_unit)

law = Eq(radius_of_orbit, first_constant + second_constant * 2**number_of_planet)


def print_law() -> str:
    return print_expression(law)


@validate_input(number_of_planet_=number_of_planet)
@validate_output(radius_of_orbit)
def calculate_radius_of_orbit(number_of_planet_: float) -> Quantity:
    if number_of_planet_ < -1:
        raise ValueError(
            "The planet number must be greater than -1 or equal."
        )

    result_expr = solve(law, radius_of_orbit, dict=True)[0][radius_of_orbit]
    result_expr = result_expr.subs({
        number_of_planet: number_of_planet_,
    })
    return Quantity(result_expr)

from sympy import (Eq, solve, sqrt)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## The third cosmic velocity is the minimum velocity that must be given to a body located near the Earth's surface so
## that it can overcome the gravitational attraction of the Earth and the Sun and leave the Solar System.

## Law is: v3 = sqrt((sqrt(2) - 1)^2 * v^2 + v2^2), where
## v3 - third cosmic velocity for a given planet,
## v - orbital velocity of the planet relative to the Sun,
## v2 - second cosmic velocity for a given planet.

third_velocity = Symbol("third_velocity", units.velocity)

orbital_velocity = Symbol("orbital_velocity", units.velocity)
second_velocity = Symbol("second_velocity", units.velocity)

law = Eq(third_velocity, sqrt(((sqrt(2) - 1)**2) * orbital_velocity**2 + second_velocity**2))


def print_law() -> str:
    return print_expression(law)


@validate_input(orbital_velocity_=orbital_velocity, second_velocity_=second_velocity)
@validate_output(third_velocity)
def calculate_third_velocity(orbital_velocity_: Quantity, second_velocity_: Quantity) -> Quantity:
    result_expr = solve(law, third_velocity, dict=True)[0][third_velocity]
    result_expr = result_expr.subs({
        orbital_velocity: orbital_velocity_,
        second_velocity: second_velocity_,
    })
    return Quantity(result_expr)

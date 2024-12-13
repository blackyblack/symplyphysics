from sympy import (Eq, solve, sqrt)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
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

# Links: Wikipedia <https://ru.wikipedia.org/wiki/%D0%A2%D1%80%D0%B5%D1%82%D1%8C%D1%8F_%D0%BA%D0%BE%D1%81%D0%BC%D0%B8%D1%87%D0%B5%D1%81%D0%BA%D0%B0%D1%8F_%D1%81%D0%BA%D0%BE%D1%80%D0%BE%D1%81%D1%82%D1%8C#cite_note-3>
# TODO: find English link

third_velocity = Symbol("third_velocity", units.velocity)

orbital_velocity = Symbol("orbital_velocity", units.velocity)
second_velocity = Symbol("second_velocity", units.velocity)

law = Eq(third_velocity, sqrt(((sqrt(2) - 1)**2) * orbital_velocity**2 + second_velocity**2))


@validate_input(orbital_velocity_=orbital_velocity, second_velocity_=second_velocity)
@validate_output(third_velocity)
def calculate_third_velocity(orbital_velocity_: Quantity, second_velocity_: Quantity) -> Quantity:
    result_expr = solve(law, third_velocity, dict=True)[0][third_velocity]
    result_expr = result_expr.subs({
        orbital_velocity: orbital_velocity_,
        second_velocity: second_velocity_,
    })
    return Quantity(result_expr)

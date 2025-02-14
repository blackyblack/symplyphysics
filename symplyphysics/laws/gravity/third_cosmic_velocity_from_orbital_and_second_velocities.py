"""
Third cosmic speed from orbital and second cosmic speed
=======================================================

The third cosmic velocity is the minimum velocity that must be given to a body located near the Earth's surface so
that it can overcome the gravitational attraction of the Earth and the Sun and leave the Solar System.

**Links:**

#. `Wikipedia <https://ru.wikipedia.org/wiki/%D0%A2%D1%80%D0%B5%D1%82%D1%8C%D1%8F_%D0%BA%D0%BE%D1%81%D0%BC%D0%B8%D1%87%D0%B5%D1%81%D0%BA%D0%B0%D1%8F_%D1%81%D0%BA%D0%BE%D1%80%D0%BE%D1%81%D1%82%D1%8C#cite_note-3>`__.

..
    TODO: find English link
    TODO: move to `gravity`?
"""

from sympy import (Eq, solve, sqrt)
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

third_cosmic_speed = clone_as_symbol(symbols.speed, subscript="3")
"""
Third cosmic :symbols:`speed` of the orbiting body.
"""

orbital_speed = symbols.speed
"""
Orbital :symbols:`speed` of the body around the attracting mass.
"""

second_cosmic_speed = clone_as_symbol(symbols.speed, subscript="2")
"""
Second cosmic :symbols:`speed` of the orbiting body.
"""

law = Eq(third_cosmic_speed, sqrt(((sqrt(2) - 1)**2) * orbital_speed**2 + second_cosmic_speed**2))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(orbital_velocity_=orbital_speed, second_velocity_=second_cosmic_speed)
@validate_output(third_cosmic_speed)
def calculate_third_velocity(orbital_velocity_: Quantity, second_velocity_: Quantity) -> Quantity:
    result_expr = solve(law, third_cosmic_speed, dict=True)[0][third_cosmic_speed]
    result_expr = result_expr.subs({
        orbital_speed: orbital_velocity_,
        second_cosmic_speed: second_velocity_,
    })
    return Quantity(result_expr)

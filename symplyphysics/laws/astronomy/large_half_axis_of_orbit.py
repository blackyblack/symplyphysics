"""
Semimajor axis of orbit via mass and speed
==========================================

Let the body move in an elliptical orbit. Then its semimajor axis depends on the mass of the body around which it
rotates and the orbital velocity.

**Notation:**

#. :quantity_notation:`gravitational_constant`.

..
    TODO: find link
    TODO: move to `gravity`?
    TODO: rename file
"""

from sympy import Eq, solve
from symplyphysics import (
    symbols,
    Quantity,
    validate_input,
    validate_output,
    quantities,
)

semimajor_axis = symbols.semimajor_axis
"""
:symbols:`semimajor_axis` of the orbit.
"""

orbital_speed = symbols.speed
"""
:symbols:`speed` of the attracted body at which it orbits the attracting body.
"""

attracting_mass = symbols.mass
"""
:symbols:`mass` of the attracting body.
"""

law = Eq(semimajor_axis, quantities.gravitational_constant * attracting_mass / orbital_speed**2)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(orbital_speed_=orbital_speed, planet_mass_=attracting_mass)
@validate_output(semimajor_axis)
def calculate_large_half_axis_length(orbital_speed_: Quantity, planet_mass_: Quantity) -> Quantity:
    result_expr = solve(law, semimajor_axis, dict=True)[0][semimajor_axis]
    result_expr = result_expr.subs({
        orbital_speed: orbital_speed_,
        attracting_mass: planet_mass_,
    })
    return Quantity(result_expr)

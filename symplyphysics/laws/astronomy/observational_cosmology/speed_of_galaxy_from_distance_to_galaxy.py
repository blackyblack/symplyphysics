"""
Speed of galaxy from distance to galaxy
=======================================

Hubble's Law is a cosmological law describing the expansion of the universe, according
to which the recessional speed of a galaxy is proportional to the distance to it from
the observer.

**Notation:**

#. :quantity_notation:`hubble_constant`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Hubble%27s_law#>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    quantities,
)

recessional_speed = symbols.speed
"""
Recessional :symbols:`speed` of the galaxy, i.e. the speed at which it moves away from the observer.
"""

distance = symbols.euclidean_distance
"""
Proper :symbols:`euclidean_distance` between the observer and the galaxy.
"""

law = Eq(recessional_speed, quantities.hubble_constant * distance)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(distance_to_galaxy_=distance)
@validate_output(recessional_speed)
def calculate_speed(distance_to_galaxy_: Quantity) -> Quantity:
    result_expr = solve(law, recessional_speed, dict=True)[0][recessional_speed]
    result_expr = result_expr.subs({
        distance: distance_to_galaxy_,
    })
    return Quantity(result_expr)

"""
Spacetime interval via time and distance
========================================

The spacetime interval is a combination of distance and time that is invariant under Lorentz tranformations.
It has the property of being invariant in the sense that it has the same value for all observers in any
inertial reference frame.

**Notation:**

#. :quantity_notation:`speed_of_light`.

**Notes:**

#. If :math:`\\Delta s**2 > 0`, the spacetime interval is said to be *timelike*. Events with a timelike
   separation can be causally connected, i.e. one can find an inertial reference frame in which both
   events happen at the same place at different times.
#. If :math:`\\Delta s**2 = 0`, the spacetime interval is said to be *lightlike*. Events with a lightlike
   separation are exactly far enough from each other that light could be present at both events, and they
   are causally connected.
#. If :math:`\\Delta s**2 < 0`, the spacetime interval is said to be *spacelike*. Events with a spacelike
   separation are causally disconnected, i.e. one can find an inertial reference frame in which both
   events happen at the same time in different positions.

**Conditions:**

#. The spacetime in which the two events occur is flat (special relativity case).

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Spacetime#Spacetime_interval>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    quantities,
)

spacetime_interval = symbols.spacetime_interval
"""
:symbols:`spacetime_interval` between the two events.
"""

temporal_distance = symbols.time
"""
:symbols:`time` separation between the two events.
"""

spatial_distance = symbols.euclidean_distance
"""
:symbols:`euclidean_distance` between the two events.
"""

law = Eq(spacetime_interval**2,
    (quantities.speed_of_light * temporal_distance)**2 - spatial_distance**2)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    temporal_distance_=temporal_distance,
    spatial_distance_=spatial_distance,
)
@validate_output(spacetime_interval)
def calculate_spacetime_interval(
    temporal_distance_: Quantity,
    spatial_distance_: Quantity,
) -> Quantity:
    expr = solve(law, spacetime_interval)[1]
    result = expr.subs({
        temporal_distance: temporal_distance_,
        spatial_distance: spatial_distance_,
    })
    return Quantity(result)

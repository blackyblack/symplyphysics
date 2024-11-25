"""
Easterly deviation from plumbline of falling bodies
===================================================

Suppose a body is falling freely in Earth's gravity field with its initial velocity being zero.
Then the effect of the Coriolis force on the falling body can be found in the fact that it deflects
from plumbline in the easterly and southerly (equatorial) directions.

**Conditions:**

#. The vector of free fall acceleration is constant.
"""

from sympy import Eq, cos, pi
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)
from symplyphysics.core.symbols.quantities import scale_factor

easterly_deviation_from_plumbline = clone_as_symbol(
    symbols.euclidean_distance,
    display_symbol="s_east",
    display_latex="s_\\text{east}",
)
"""
Easterly deviation of falling body from plumbline due to Earth's rotation.
See :symbols:`euclidean_distance`
"""

fall_time = symbols.time
"""
:symbols:`time` of the body's fall.
"""

rotation_period = symbols.period
"""
:symbols:`period` of the Earth's rotation.
"""

initial_elevation = symbols.height
"""
Initial elevation (:symbols:`height`) of the body from the Earth's surface.
"""

latitude = symbols.latitude
"""
:symbols:`latitude` of the location of the body.
"""

law = Eq(
    easterly_deviation_from_plumbline,
    (4 * pi / 3) * (fall_time / rotation_period) * initial_elevation * cos(latitude),
)
"""
:laws:symbol::

:laws:latex::
"""

# TODO: derive from the solution of the relative motion equation in Earth's gravitational field.


@validate_input(
    fall_time_=fall_time,
    rotation_period_=rotation_period,
    initial_elevation_=initial_elevation,
    latitude_=latitude,
)
@validate_output(easterly_deviation_from_plumbline)
def calculate_easterly_deviation_from_plumbline(
    fall_time_: Quantity,
    rotation_period_: Quantity,
    initial_elevation_: Quantity,
    latitude_: Quantity | float,
) -> Quantity:
    result = law.rhs.subs({
        fall_time: fall_time_,
        rotation_period: rotation_period_,
        initial_elevation: initial_elevation_,
        latitude: scale_factor(latitude_),
    })
    return Quantity(result)

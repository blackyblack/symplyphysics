"""
Southerly deviation from plumbline of falling bodies
====================================================

Suppose a body is falling freely in Earth's gravity field with its initial velocity being zero.
Then the effect of the Coriolis force on the falling body can be found in the fact that it deflects
from plumbline in the easterly and southerly (equatorial) directions.

**Notes:**

#. The southerly deviation is extremely small and almost unobservable due to the :math:`t / T` factor.

#. Also see :ref:`Easterly deviation from plumbline of falling bodies`.

**Conditions:**

#. The vector of free fall acceleration is considered constant.

**Links:**

#. Sivukhin D.V. (1979), __Obshchiy kurs fiziki__ [General course of Physics], vol. 1, p. 355, (67.11).
"""

from sympy import Eq, sin, pi
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)
from symplyphysics.core.symbols.quantities import scale_factor

southerly_deviation_from_plumbline = clone_as_symbol(
    symbols.euclidean_distance,
    display_symbol="s_south",
    display_latex="s_\\text{south}",
)
"""
Southerly (equatorial) deviation of falling body from plumbline due to Earth's rotation.
See :symbols:`euclidean_distance`.
"""

fall_time = symbols.time
"""
:symbols:`time` elapsed during the body's fall.
"""

rotation_period = symbols.period
"""
:symbols:`period` of Earth's rotation.
"""

easterly_deviation_from_plumbline = clone_as_symbol(
    symbols.euclidean_distance,
    display_symbol="s_east",
    display_latex="s_\\text{east}",
)
"""
Easterly deviation of falling body from plumbline. See :symbols:`euclidean_distance`.
"""

latitude = symbols.latitude
"""
:symbols:`latitude` of the place the body is located in.
"""

law = Eq(
    southerly_deviation_from_plumbline,
    pi * (fall_time / rotation_period) * easterly_deviation_from_plumbline * sin(latitude),
)
"""
:laws:symbol::

:laws:latex::
"""

# TODO: derive from the solution of the relative motion equation in Earth's gravitational field.


@validate_input(
    fall_time_=fall_time,
    rotation_period_=rotation_period,
    easterly_deviation_from_plumbline_=easterly_deviation_from_plumbline,
    latitude_=latitude,
)
@validate_output(easterly_deviation_from_plumbline)
def calculate_southerly_deviation_from_plumbline(
    fall_time_: Quantity,
    rotation_period_: Quantity,
    easterly_deviation_from_plumbline_: Quantity,
    latitude_: Quantity | float,
) -> Quantity:
    result = law.rhs.subs({
        fall_time: fall_time_,
        rotation_period: rotation_period_,
        easterly_deviation_from_plumbline: easterly_deviation_from_plumbline_,
        latitude: scale_factor(latitude_),
    })
    return Quantity(result)

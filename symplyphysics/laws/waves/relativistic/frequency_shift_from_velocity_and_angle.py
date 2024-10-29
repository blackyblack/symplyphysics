"""
Frequency shift from speed and angle
====================================

See :doc:`laws.waves.relativistic.longitudinal_frequency_shift_from_velocity`.

**Notes:**

#. It is not trivial to substitute moving source and idle observer with moving observer and idle source. Law
   depends on the current frame (observer or source are at rest at this frame), angle is detected at the point
   of emission or at the point of reception (see relativistic angle aberration).

**Conditions:**

#. Angle is measured at the moment of emission with respect to the observer frame.
#. Motion is in 2D space.
"""

from sympy import Eq, cos, solve, sqrt
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    clone_as_symbol,
    symbols,
)
from symplyphysics.quantities import speed_of_light
from symplyphysics.core.symbols.quantities import scale_factor

observer_frequency = clone_as_symbol(
    symbols.temporal_frequency,
    display_symbol="f_o",
    display_latex="f_\\text{o}",
)
"""
Observed :symbols:`temporal_frequency` of the wave.
"""

source_frequency = clone_as_symbol(
    symbols.temporal_frequency,
    display_symbol="f_s",
    display_latex="f_\\text{s}",
)
"""
Source :symbols:`temporal_frequency` of the wave.
"""

relative_speed = symbols.speed
"""
Relative :symbols:`speed` between the source and the observer.
"""

source_angle = symbols.angle
"""
:symbols:`angle` between the signal vector (directed from the source to the observer) and the source's velocity
vector.
"""

law = Eq(
    observer_frequency,
    source_frequency * sqrt(speed_of_light**2 - relative_speed**2) /
    (speed_of_light - relative_speed * cos(source_angle)))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(real_frequency_=source_frequency,
    relative_speed_=relative_speed,
    source_angle_=source_angle)
@validate_output(observer_frequency)
def calculate_observed_frequency(real_frequency_: Quantity, relative_speed_: Quantity,
    source_angle_: float | Quantity) -> Quantity:
    #HACK: sympy angles are always in radians
    source_angle_radians = scale_factor(source_angle_)
    result_expr = solve(law, observer_frequency, dict=True)[0][observer_frequency]
    frequency_applied = result_expr.subs({
        source_frequency: real_frequency_,
        relative_speed: relative_speed_,
        source_angle: source_angle_radians
    })
    return Quantity(frequency_applied)

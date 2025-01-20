"""
Longitudinal frequency shift from speeds
========================================

General relativistic Doppler effect that is classical Doppler effect with relativistic coefficient.
This law is not used for actual calculations because relativistic effects are not visible for
acoustic waves. And for electromagnetic waves it is hard to define velocity relative to medium.
This law is used to show the connection between classical and relativistic Doppler laws.

**Notes:**

#. The speed of light is used in this law, so replace it with the correct value of the speed of
   light in the medium.
#. When wave velocity is getting close to speed of light, we are no longer having observer and source
   velocities, but relativistic relative velocity.

**Conditions:**

#. Motion is in 1D space.
#. Coordinate system is at rest with respect to the medium of the wave, or any coordinate system
   for electromagnetic wave in vacuum.
#. Medium is not fixed in space.
#. Source and observer velocities are no greater than wave velocity.
#. The source and observer velocities are collinear.

..
    TODO rename file
    TODO add link
"""

from sympy import Eq, solve, sqrt
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)
from symplyphysics.quantities import speed_of_light

observer_frequency = clone_as_symbol(
    symbols.temporal_frequency,
    display_symbol="f_o",
    display_latex="f_\\text{o}",
)
"""
Observer :symbols:`temporal_frequency`.
"""

source_frequency = clone_as_symbol(
    symbols.temporal_frequency,
    display_symbol="f_s",
    display_latex="f_\\text{s}",
)
"""
Source :symbols:`temporal_frequency`.
"""

source_speed = clone_as_symbol(
    symbols.speed,
    display_symbol="v_s",
    display_latex="v_\\text{s}",
)
"""
Source :symbols:`speed`, positive when moving away from observer and negative otherwise.
"""

observer_speed = clone_as_symbol(
    symbols.speed,
    display_symbol="v_o",
    display_latex="v_\\text{o}",
)
"""
Observer :symbols:`speed`, positive when moving away from source and negative otherwise.
"""

wave_speed = symbols.phase_speed
"""
:symbols:`phase_speed` of the wave.
"""

law = Eq(
    observer_frequency,
    source_frequency * (1 - observer_speed / wave_speed) / (1 + source_speed / wave_speed) * sqrt(
    (1 - (source_speed / speed_of_light)**2) / (1 - (observer_speed / speed_of_light)**2)))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(real_frequency_=source_frequency,
    wave_velocity_=wave_speed,
    source_velocity_=source_speed,
    observer_velocity_=observer_speed)
@validate_output(observer_frequency)
def calculate_observed_frequency(real_frequency_: Quantity, wave_velocity_: Quantity,
    source_velocity_: Quantity, observer_velocity_: Quantity) -> Quantity:

    result_expr = solve(law, observer_frequency, dict=True)[0][observer_frequency]
    frequency_applied = result_expr.subs({
        source_frequency: real_frequency_,
        wave_speed: wave_velocity_,
        source_speed: source_velocity_,
        observer_speed: observer_velocity_
    })
    return Quantity(frequency_applied)

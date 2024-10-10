"""
Frequency shift from speed in collinear motion
==============================================

The *Doppler effect* or *Doppler shift* is the apparent change in frequency of a wave in
relation to an observer moving relative to the wave source.

**Notes:**

#. When the source is moving with speed equal to wave speed, the observed frequency goes
   to infinity. This effect is known as the *sonic boom*.

**Conditions:**

#. The source and observer speeds are less or equal to the wave speed. Otherwise emitted
   waves are left behind the source or never reach the observer.
#. The source and observer are moving directly towards or away from each other (collinear
   motion).
"""

from sympy import Eq, pi, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
    quantities,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.waves import (
    frequency_shift_from_speed_in_arbitrary_motion as classical_doppler_with_angle,
)
from symplyphysics.laws.waves.relativistic import (
    longitudinal_frequency_shift_from_absolute_velocities as general_doppler_law,
)

observed_frequency = clone_as_symbol(
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
Wave :symbols:`temporal_frequency` of the source.
"""

wave_speed = symbols.phase_speed    
"""
:symbols:`phase_speed` of the wave.
"""

source_speed = clone_as_symbol(
    symbols.speed,
    display_symbol="v_s",
    display_latex="v_\\text{s}",
)
"""
:symbols:`speed` of the wave source, positive when moving away from the observer and negative otherwise.
"""

observer_speed = clone_as_symbol(
    symbols.speed,
    display_symbol="v_o",
    display_latex="v_\\text{o}",
)
"""
:symbols:`speed` of the observer, positive when moving away from the source and negative otherwise.
"""

law = Eq(observed_frequency,
    source_frequency * (wave_speed - observer_speed) / (wave_speed + source_speed))
"""
:laws:symbol::

:laws:latex::
"""

# Confirm that classical Doppler effect is a special case of relativistic Doppler effect

# This is a general formula for Doppler effect that has both classical and relativistic parts of equation.

# Relativistic part is zero for velocities much less than speed of light
_classical_law = general_doppler_law.law.subs({
    general_doppler_law.source_frequency: source_frequency,
    general_doppler_law.wave_speed: wave_speed,
    general_doppler_law.source_speed: source_speed,
    general_doppler_law.observer_speed: observer_speed,
    (general_doppler_law.source_speed / quantities.speed_of_light)**2: 0,
    (general_doppler_law.observer_speed / quantities.speed_of_light)**2: 0
})
assert expr_equals(_classical_law.rhs, law.rhs)

# Confirm that Doppler effect for collinear movement is a subset of Doppler effect with angles

## Classical Doppler effect angles are calculated with respect to the signal vector, directed
## from source to observer. Hence source moving directly towards observer has 0 angle. Therefore
## its cosine is positive, when moving towards observer.
## This law has reverse notation - source velocity is positive when moving away from the observer.
## Therefore we should use opposite direction for source - set pi as source angle.
_observed_frequency_zero_angles = classical_doppler_with_angle.law.subs({
    classical_doppler_with_angle.observer_angle: 0,
    classical_doppler_with_angle.source_angle: pi,
    classical_doppler_with_angle.observer_speed: observer_speed,
    classical_doppler_with_angle.source_speed: source_speed,
    classical_doppler_with_angle.source_frequency: source_frequency,
    classical_doppler_with_angle.wave_speed: wave_speed
}).rhs

assert expr_equals(_observed_frequency_zero_angles, law.rhs)


@validate_input(real_frequency_=source_frequency,
    wave_velocity_=wave_speed,
    source_velocity_=source_speed,
    observer_velocity_=observer_speed)
@validate_output(observed_frequency)
def calculate_observed_frequency(real_frequency_: Quantity, wave_velocity_: Quantity,
    source_velocity_: Quantity, observer_velocity_: Quantity) -> Quantity:
    result_expr = solve(law, observed_frequency, dict=True)[0][observed_frequency]
    frequency_applied = result_expr.subs({
        source_frequency: real_frequency_,
        wave_speed: wave_velocity_,
        source_speed: source_velocity_,
        observer_speed: observer_velocity_
    })
    return Quantity(frequency_applied)

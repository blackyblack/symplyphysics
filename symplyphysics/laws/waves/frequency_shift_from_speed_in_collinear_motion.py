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

from sympy import (Eq, pi, solve)
from sympy.physics.units import speed_of_light
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.waves import frequency_shift_from_velocity_and_angle as classical_doppler_with_angle
from symplyphysics.laws.waves.relativistic import longitudinal_frequency_shift_from_absolute_velocities as general_doppler_law

observed_frequency = Symbol("observed_frequency", units.frequency)
r"""
Observed frequency of the wave.

Symbol:
    :code:`f_o`

Latex:
    :math:`f_\text{o}`
"""

source_frequency = Symbol("source_frequency", units.frequency)
r"""
Wave frequency of the source.

Symbol:
    :code:`f_s`

Latex:
    :math:`f_\text{s}`
"""

wave_speed = Symbol("wave_speed", units.velocity)
"""
Phase speed of the wave.

Symbol:
    :code:`v`
"""

source_speed = Symbol("source_speed", units.velocity)
r"""
Speed of the wave source, positive when moving away from the observer and negative otherwise.

Symbol:
    :code:`v_s`

Latex:
    :math:`v_\text{s}`
"""

observer_speed = Symbol("observer_speed", units.velocity)
r"""
Speed of the observer, positive when moving away from the source and negative otherwise.

Symbol:
    :code:`v_o`

Latex:
    :math:`v_\text{o}`
"""

law = Eq(observed_frequency,
    source_frequency * (wave_speed - observer_speed) / (wave_speed + source_speed))
r"""
:code:`f_o = f_s * (v - v_o) / (v + v_s)`

Latex:
    .. math::
        f_\text{o} = f_\text{s} \frac{v - v_\text{o}}{v + v_\text{s}}
"""

# Confirm that classical Doppler effect is a special case of relativistic Doppler effect

# This is a general formula for Doppler effect that has both classical and relativistic parts of equation.

# Relativistic part is zero for velocities much less than speed of light
classical_law = general_doppler_law.law.subs({
    general_doppler_law.real_frequency: source_frequency,
    general_doppler_law.wave_velocity: wave_speed,
    general_doppler_law.source_velocity: source_speed,
    general_doppler_law.observer_velocity: observer_speed,
    (general_doppler_law.source_velocity / speed_of_light)**2: 0,
    (general_doppler_law.observer_velocity / speed_of_light)**2: 0
})
assert expr_equals(classical_law.rhs, law.rhs)

# Confirm that Doppler effect for collinear movement is a subset of Doppler effect with angles

## Classical Doppler effect angles are calculated with respect to the signal vector, directed
## from source to observer. Hence source moving directly towards observer has 0 angle. Therefore
## its cosine is positive, when moving towards observer.
## This law has reverse notation - source velocity is positive when moving away from the observer.
## Therefore we should use opposite direction for source - set pi as source angle.
observed_frequency_zero_angles = classical_doppler_with_angle.law.subs({
    classical_doppler_with_angle.observer_angle: 0,
    classical_doppler_with_angle.source_angle: pi,
    classical_doppler_with_angle.observer_speed: observer_speed,
    classical_doppler_with_angle.source_speed: source_speed,
    classical_doppler_with_angle.real_frequency: source_frequency,
    classical_doppler_with_angle.wave_velocity: wave_speed
}).rhs

assert expr_equals(observed_frequency_zero_angles, law.rhs)


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

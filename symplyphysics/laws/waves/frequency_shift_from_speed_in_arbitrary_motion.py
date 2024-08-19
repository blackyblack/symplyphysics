"""
Frequency shift from speed in arbitrary motion
==============================================

The *Doppler effect* or *Doppler shift* is the apparent change in frequency of a wave in
relation to an observer moving relative to the wave source.

Also see :doc:`laws.waves.frequency_shift_from_speed_in_collinear_motion`.

**Conditions:**

#. The source and observer speeds are less or equal to the wave speed. Otherwise emitted
   waves are left behind the source or never reach the observer.
#. The speeds are much less than the speed of light, i.e. this law is non-relativistic.
"""

from sympy import (Eq, cos, solve)
from symplyphysics import (units, angle_type, Quantity, Symbol, validate_input,
    validate_output)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.symbols.quantities import scale_factor
from symplyphysics.laws.waves import wavelength_from_phase_speed_and_period as period_law
from symplyphysics.definitions import temporal_frequency_from_period as frequency_def
from symplyphysics.laws.geometry import planar_projection_is_cosine as projector
from symplyphysics.laws.kinematics import distance_from_constant_velocity as distance_law

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
Magnitude of the velocity vector of the source.

Symbol:
    :code:`v_s`

Latex:
    :math:`v_\text{s}`
"""

observer_speed = Symbol("observer_speed", units.velocity)
r"""
Magnitude of the velocity vector of the observer.

Symbol:
    :code:`v_o`

Latex:
    :math:`v_\text{o}`
"""

source_angle = Symbol("source_angle", angle_type)
r"""
Angle between the wave velocity and the source velocity.

Symbol:
    :code:`theta_s`

Latex:
    :math:`\theta_\text{s}`
"""

observer_angle = Symbol("observer_angle", angle_type)
r"""
Angle between the wave velocity and the observer velocity.

Symbol:
    :code:`theta_o`

Latex:
    :math:`\theta_\text{o}`
"""

law = Eq(
    observed_frequency,
    source_frequency * (wave_speed - observer_speed * cos(observer_angle)) /
    (wave_speed - source_speed * cos(source_angle)))
r"""
:code:`f_o = f_s * (v - v_o * cos(theta_o)) / (v - v_s * cos(theta_s))`

Latex:
    .. math::
        f_\text{o} = f_\text{s} \frac{v - v_\text{o} \cos{\theta_\text{o}}}{v - v_\text{s} \cos{\theta_\text{s}}}
"""

# Derive the same law from frequency, wavelength laws, and assumption that moving source or
# observer affects wavelength.

## Start with idle observer and moving source

_period_from_frequency = solve(frequency_def.law, frequency_def.period,
    dict=True)[0][frequency_def.period]
_wave_period = _period_from_frequency.subs(frequency_def.temporal_frequency, source_frequency)

_wavelength_from_period = solve(period_law.law, period_law.wavelength,
    dict=True)[0][period_law.wavelength]
_wavelength = _wavelength_from_period.subs({
    period_law.period: _wave_period,
    period_law.phase_velocity: wave_speed
})

## While wave travels (_wave_period * wave_speed) distance, moving source travels (_wave_period * source_velocity)
## distance.
## We are only interested in the wavelength on the wave signal vector (shortest path from source to observer), as
## it is what we measure on observer. Therefore we take 'source_speed' projection on the signal vector.

_source_speed_projection_on_signal = projector.law.subs({
    projector.vector_angle: source_angle,
    projector.vector_length: source_speed
}).rhs

## Assume constant velocity during '_wave_period'
_moving_source_distance_for_period = distance_law.law.subs({
    distance_law.initial_position: 0,
    distance_law.movement_time: _wave_period,
    distance_law.constant_velocity: _source_speed_projection_on_signal,
}).rhs

## Assuming signal vector pointing from source to observer, positive projection should decrease _wavelength.
_wavelength_observed = _wavelength - _moving_source_distance_for_period

_period_from_wavelength = solve(period_law.law, period_law.period,
    dict=True)[0][period_law.period]
_observed_wave_period = _period_from_wavelength.subs({
    period_law.wavelength: _wavelength_observed,
    period_law.phase_velocity: wave_speed
})

## Confirm that derived law is the same as expected for idle observer

_frequency_from_period = solve(frequency_def.law, frequency_def.temporal_frequency,
    dict=True)[0][frequency_def.temporal_frequency]
_frequency_observed = _frequency_from_period.subs(frequency_def.period, _observed_wave_period)
assert expr_equals(_frequency_observed, law.rhs.subs(observer_speed, 0))

## Now apply movement of the observer

observer_speed_projection_on_signal = projector.law.subs({
    projector.vector_angle: observer_angle,
    projector.vector_length: observer_speed
}).rhs

# NOTE: Relativistic velocity addition should be applied when wave speed is close to speed of light

## Assuming signal vector pointing from source to observer, positive projection should decrease relative wave
## velocity from observer point of view, according to Galilean velocity addition formula.
_relative_wave_speed = wave_speed - observer_speed_projection_on_signal

_period_relative_source = _period_from_wavelength.subs({
    period_law.wavelength: _wavelength_observed,
    period_law.phase_velocity: _relative_wave_speed
})

_frequency_relative_observer = _frequency_from_period.subs(frequency_def.period,
    _period_relative_source)

## Confirm that derived law is the same as expected
assert expr_equals(_frequency_relative_observer, law.rhs)


@validate_input(real_frequency_=source_frequency,
    wave_velocity_=wave_speed,
    source_speed_angle=(source_speed, source_angle),
    observer_speed_angle=(observer_speed, observer_angle))
@validate_output(observed_frequency)
def calculate_observed_frequency(
        real_frequency_: Quantity, wave_velocity_: Quantity, source_speed_angle: tuple[Quantity,
    float | Quantity], observer_speed_angle: tuple[Quantity, float | Quantity]) -> Quantity:
    (observer_speed_, observer_angle_) = observer_speed_angle
    (source_speed_, source_angle_) = source_speed_angle
    #HACK: sympy angles are always in radians
    observer_angle_radians = scale_factor(observer_angle_)
    source_angle_radians = scale_factor(source_angle_)
    result_expr = solve(law, observed_frequency, dict=True)[0][observed_frequency]
    frequency_applied = result_expr.subs({
        source_frequency: real_frequency_,
        wave_speed: wave_velocity_,
        source_speed: source_speed_,
        observer_speed: observer_speed_,
        source_angle: source_angle_radians,
        observer_angle: observer_angle_radians
    })
    return Quantity(frequency_applied)

from sympy import (Eq, cos, solve)
from symplyphysics import (units, angle_type, Quantity, Symbol, print_expression, validate_input,
    validate_output)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.waves import wavelength_from_wave_speed_and_period as period_law
from symplyphysics.laws.kinematic import temporal_frequency_from_period as frequency_def
from symplyphysics.laws.kinematic import planar_projection_is_cosine as projector
from symplyphysics.laws.kinematic import distance_from_constant_velocity as distance_law

# Description
## See [doppler effect](./frequency_shift_from_velocity.py) description. When objects are not moving collinear, one
## should add angles to the formula.

# Law: fo = fs * (v - vo * cos(pho))/(v - vs * cos(phs)), where
## fo is observed frequency,
## fs is source wave frequency,
## v is wave velocity in this media,
## vo is observer speed relative to the medium (magnitude of the velocity vector),
## vs is source speed relative to the medium (magnitude of the velocity vector),
## pho is angle between vector pointing from source of the wave to the observer (signal vector) and observer velocity vector,
## phs is angle between signal vector and source velocity vector.

# Conditions:
## - Source and observer velocities are less or equal than wave velocity. Otherwise emitted waves are left behind the source or never
## reach the observer.
## - Motion is in 2-D space.
## - Source speed is constant during one period of the wave.
## - Measurements are made in the medium reference frame (medium is not moving).

# Note:
## Signal vector connects the source at the moment of wave emission and observer at the moment
## of wave reception. It may be quite hard to estimate the point where the observer will be at the
## moment of observation. Therefore this law is usually applied when objects are moving slowly, or when
## one of the objects is idle.

# Note:
## When there is no medium to choose as reference point, there are no other waves but electromagnetic
## that are known to always travel at the speed of light.
## Use [relativistic Doppler law](..\relativistic\waves\frequency_shift_from_velocity_and_angle.py)
## if no medium is present.

observed_frequency = Symbol("observed_frequency", units.frequency)
real_frequency = Symbol("real_frequency", units.frequency)
wave_velocity = Symbol("wave_velocity", units.velocity)
source_speed = Symbol("source_speed", units.velocity)
observer_speed = Symbol("observer_speed", units.velocity)
source_angle = Symbol("source_angle", angle_type)
observer_angle = Symbol("observer_angle", angle_type)

law = Eq(
    observed_frequency,
    real_frequency * (wave_velocity - observer_speed * cos(observer_angle)) /
    (wave_velocity - source_speed * cos(source_angle)))

# Derive the same law from frequency, wavelength laws, and assumption that moving source or
# observer affects wavelength.

## Start with idle observer and moving source

period_from_frequency = solve(frequency_def.law, frequency_def.period,
    dict=True)[0][frequency_def.period]
wave_period = period_from_frequency.subs(frequency_def.temporal_frequency, real_frequency)

wavelength_from_period = solve(period_law.law, period_law.wavelength,
    dict=True)[0][period_law.wavelength]
wavelength = wavelength_from_period.subs({
    period_law.oscillation_period: wave_period,
    period_law.propagation_speed: wave_velocity
})

## While wave travels (wave_period * wave_velocity) distance, moving source travels (wave_period * source_velocity)
## distance.
## We are only interested in the wavelength on the wave signal vector (shortest path from source to observer), as
## it is what we measure on observer. Therefore we take 'source_speed' projection on the signal vector.

source_speed_projection_on_signal = projector.law.subs({
    projector.vector_angle: source_angle,
    projector.vector_length: source_speed
}).rhs

## Assume constant velocity during 'wave_period'
moving_source_distance_for_period = distance_law.law.subs({
    distance_law.initial_position: 0,
    distance_law.movement_time: wave_period,
    distance_law.constant_velocity: source_speed_projection_on_signal,
}).rhs

## Assuming signal vector pointing from source to observer, positive projection should decrease wavelength.
wavelength_observed = wavelength - moving_source_distance_for_period

period_from_wavelength = solve(period_law.law, period_law.oscillation_period,
    dict=True)[0][period_law.oscillation_period]
observed_wave_period = period_from_wavelength.subs({
    period_law.wavelength: wavelength_observed,
    period_law.propagation_speed: wave_velocity
})

## Confirm that derived law is the same as expected for idle observer

frequency_from_period = solve(frequency_def.law, frequency_def.temporal_frequency,
    dict=True)[0][frequency_def.temporal_frequency]
frequency_observed = frequency_from_period.subs(frequency_def.period, observed_wave_period)
assert expr_equals(frequency_observed, law.rhs.subs(observer_speed, 0))

## Now apply movement of the observer

observer_speed_projection_on_signal = projector.law.subs({
    projector.vector_angle: observer_angle,
    projector.vector_length: observer_speed
}).rhs

# NOTE: Relativistic velocity addition should be applied when wave speed is close to speed of light

## Assuming signal vector pointing from source to observer, positive projection should decrease relative wave
## velocity from observer point of view, according to Galilean velocity addition formula.
relative_wave_speed = wave_velocity - observer_speed_projection_on_signal

period_relative_source = period_from_wavelength.subs({
    period_law.wavelength: wavelength_observed,
    period_law.propagation_speed: relative_wave_speed
})

frequency_relative_observer = frequency_from_period.subs(frequency_def.period,
    period_relative_source)

## Confirm that derived law is the same as expected
assert expr_equals(frequency_relative_observer, law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(real_frequency_=real_frequency,
    wave_velocity_=wave_velocity,
    source_speed_angle=(source_speed, source_angle),
    observer_speed_angle=(observer_speed, observer_angle))
@validate_output(observed_frequency)
def calculate_observed_frequency(
        real_frequency_: Quantity, wave_velocity_: Quantity, source_speed_angle: tuple[Quantity,
    float | Quantity], observer_speed_angle: tuple[Quantity, float | Quantity]) -> Quantity:
    (observer_speed_, observer_angle_) = observer_speed_angle
    (source_speed_, source_angle_) = source_speed_angle
    #HACK: sympy angles are always in radians
    observer_angle_radians = observer_angle_.scale_factor if isinstance(observer_angle_,
        Quantity) else observer_angle_
    source_angle_radians = source_angle_.scale_factor if isinstance(source_angle_,
        Quantity) else source_angle_
    result_expr = solve(law, observed_frequency, dict=True)[0][observed_frequency]
    frequency_applied = result_expr.subs({
        real_frequency: real_frequency_,
        wave_velocity: wave_velocity_,
        source_speed: source_speed_,
        observer_speed: observer_speed_,
        source_angle: source_angle_radians,
        observer_angle: observer_angle_radians
    })
    return Quantity(frequency_applied)

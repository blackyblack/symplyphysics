import numbers
from sympy import (Eq, cos, solve)
from symplyphysics import (units, angle_type, expr_to_quantity, Quantity, Symbol, print_expression,
    validate_input_symbols, validate_output_symbol)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.waves import wavelength_from_wave_speed_and_period as period_law
from symplyphysics.laws.kinematic import temporal_frequency_from_period as frequency_def
from symplyphysics.laws.kinematic import planar_projection_is_cosine as projector

# Description
## See [doppler effect](./frequency_shift_from_velocity.py) description. When objects are not moving collinear, one
## should add angles to the formula.

# Law: fo = fs * (v - vo * cos(pho))/(v - vs * cos(phs)), where
## fo is observed frequency,
## fs is source wave frequency,
## v is wave velocity in this media,
## vo is observer speed relative to the medium (magnitude of the velocity vector),
## vs is source speed relative to the medium (magnitude of the velocity vector),
## pho is angle between signal vector and observer velocity vector,
## phs is angle between vector pointing from source of the wave to the observer (signal vector), and source velocity vector.

# Conditions:
## - Source and observer velocities are less or equal than wave velocity. Otherwise emitted waves are left behind the source or never
## reach the observer.
## - Motion is in 2-D space.
## - Angles are measured instantly for both source and observer.

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

wave_period_from_frequency_solved = solve(frequency_def.law, frequency_def.period, dict=True)[0][frequency_def.period]
wave_period = wave_period_from_frequency_solved.subs(frequency_def.temporal_frequency, real_frequency)

wavelength_solved = solve(period_law.law, period_law.wavelength, dict=True)[0][period_law.wavelength]
wavelength = wavelength_solved.subs({
    period_law.oscillation_period: wave_period,
    period_law.propagation_speed: wave_velocity
})

## While wave travels (wave_period * wave_velocity) distance, moving source travels (wave_period * source_velocity)
## distance.
## We are only interested in the wavelength on the wave signal vector (shortest path from source to observer), as
## it is what we measure. Therefore we take 'source_speed' projection on the signal vector.

# NOTE: Angles at the time of emission and at the time of observation may change. Here we assume, that wave
# propagation is so fast that angles do not change.

source_speed_projection_on_signal = projector.law.subs({
    projector.vector_angle: source_angle,
    projector.vector_length: source_speed
}).rhs

moving_source_distance_for_period = wave_period * source_speed_projection_on_signal

## Assuming signal vector pointing from source to observer, positive projection should decrease wavelength.
wavelength_observed = wavelength - moving_source_distance_for_period

wave_period_solved = solve(period_law.law, period_law.oscillation_period, dict=True)[0][period_law.oscillation_period]
observed_wave_period = wave_period_solved.subs({
    period_law.wavelength: wavelength_observed,
    period_law.propagation_speed: wave_velocity
})

frequency_solved = solve(frequency_def.law, frequency_def.temporal_frequency, dict=True)[0][frequency_def.temporal_frequency]
frequency_observed = frequency_solved.subs(frequency_def.period, observed_wave_period)

## Confirm that derived law is the same as expected for idle observer
assert expr_equals(frequency_observed, law.rhs.subs(observer_speed, 0))

observer_speed_projection_on_signal = projector.law.subs({
    projector.vector_angle: observer_angle,
    projector.vector_length: observer_speed
}).rhs

# NOTE: Relativistic velocity addition should be applied when wave speed is close to speed of light

## Assuming signal vector pointing from source to observer, positive projection should decrease relative wave
## velocity from observer point of view is lower, according to Galilean velocity addition formula.
relative_wave_speed = wave_velocity - observer_speed_projection_on_signal

period_relative_source = wave_period_solved.subs({
    period_law.wavelength: wavelength_observed,
    period_law.propagation_speed: relative_wave_speed
})

frequency_relative_observer = frequency_solved.subs(frequency_def.period, period_relative_source)

## Confirm that derived law is the same as expected
assert expr_equals(frequency_relative_observer, law.rhs)


def print() -> str:
    return print_expression(law)


@validate_input_symbols(real_frequency_=real_frequency,
    wave_velocity_=wave_velocity,
    source_speed_=source_speed,
    observer_speed_=observer_speed,
    source_angle_=source_angle,
    observer_angle_=observer_angle)
@validate_output_symbol(observed_frequency)
def calculate_observed_frequency(real_frequency_: Quantity, wave_velocity_: Quantity,
    source_speed_: Quantity, observer_speed_: Quantity, source_angle_: float | Quantity,
    observer_angle_: float | Quantity) -> Quantity:
    #HACK: sympy angles are always in radians
    source_angle_radians = source_angle_ if isinstance(source_angle_,
        numbers.Number) else source_angle_.scale_factor
    observer_angle_radians = observer_angle_ if isinstance(observer_angle_,
        numbers.Number) else observer_angle_.scale_factor
    result_expr = solve(law, observed_frequency, dict=True)[0][observed_frequency]
    frequency_applied = result_expr.subs({
        real_frequency: real_frequency_,
        wave_velocity: wave_velocity_,
        source_speed: source_speed_,
        observer_speed: observer_speed_,
        source_angle: source_angle_radians,
        observer_angle: observer_angle_radians
    })
    return expr_to_quantity(frequency_applied)

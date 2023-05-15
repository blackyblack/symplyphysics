from sympy import (Eq, solve)
from sympy.physics.units import speed_of_light
from symplyphysics import (units, expr_to_quantity, Quantity, Symbol, print_expression,
    validate_input_symbols, validate_output_symbol)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.relativistic.waves import longitudinal_frequency_shift_from_absolute_velocities as general_doppler_law

#TODO: add Doppler effect with angle

# Description
## The Doppler effect or Doppler shift is the apparent change in frequency of a wave in relation to an observer moving relative to the wave source.

# Law: fo = fs * (v - vo)/(v + vs), where
## fo is observed frequency,
## fs is source wave frequency,
## v is wave velocity in this media,
## vo is observer velocity relative to the medium (positive when moving away from source, negative when moving towards source),
## vs is source velocity relative to the medium (positive when moving away from observer, negative when moving towards observer).

# Conditions:
## - Source and observer velocities are less or equal than wave velocity. Otherwise emitted waves are left behind the source or never
## reach the observer.
## - Source and observer are moving directly towards or away from each other (collinear motion).

# Note:
## When source is moving with speed equal to wave velocity, observed frequency goes to infinity. This effect is known as sonic boom.

# Note:
## When direction of source and observer is not directed towards to/away from each other, one should introduce angle to the formula.

# Note:
## Moving source is equivalent to the moving observer in the opposite direction, if velocities are much less than wave velocity. If
## velocities are comparable to wave velocity, Doppler shift is no longer irrelative.

observed_frequency = Symbol("observed_frequency", units.frequency)
real_frequency = Symbol("real_frequency", units.frequency)
wave_velocity = Symbol("wave_velocity", units.velocity)
source_velocity = Symbol("source_velocity", units.velocity)
observer_velocity = Symbol("observer_velocity", units.velocity)

law = Eq(observed_frequency,
    real_frequency * (wave_velocity - observer_velocity) / (wave_velocity + source_velocity))

# Confirm that classical Doppler effect is a special case of relativistic Doppler effect

# This is a general formula for Doppler effect that has both classical and relativistic parts of equation.

# Relativistic part is zero for velocities much less than speed of light
classical_law = general_doppler_law.law.subs({
    general_doppler_law.real_frequency: real_frequency,
    general_doppler_law.wave_velocity: wave_velocity,
    general_doppler_law.source_velocity: source_velocity,
    general_doppler_law.observer_velocity: observer_velocity,
    (general_doppler_law.source_velocity / speed_of_light)**2: 0,
    (general_doppler_law.observer_velocity / speed_of_light)**2: 0})
assert expr_equals(classical_law.rhs, law.rhs)

#TODO: add proof for classical Doppler effect


def print() -> str:
    return print_expression(law)


@validate_input_symbols(real_frequency_=real_frequency,
    wave_velocity_=wave_velocity,
    source_velocity_=source_velocity,
    observer_velocity_=observer_velocity)
@validate_output_symbol(observed_frequency)
def calculate_observed_frequency(real_frequency_: Quantity, wave_velocity_: Quantity,
    source_velocity_: Quantity, observer_velocity_: Quantity) -> Quantity:
    result_expr = solve(law, observed_frequency, dict=True)[0][observed_frequency]
    frequency_applied = result_expr.subs({
        real_frequency: real_frequency_,
        wave_velocity: wave_velocity_,
        source_velocity: source_velocity_,
        observer_velocity: observer_velocity_
    })
    return expr_to_quantity(frequency_applied)

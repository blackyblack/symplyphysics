from sympy import (Eq, solve)
from symplyphysics import (units, expr_to_quantity, Quantity, Symbol, print_expression,
    validate_input_symbols, validate_output_symbol)

# Description
## The Doppler effect or Doppler shift is the apparent change in frequency of a wave in relation to an observer moving relative to the wave source.

# Law: fo = fs * (v + vo)/(v + vs), where
## fo is observed frequency,
## fs is source wave frequency,
## v is wave velocity in this media,
## vo is observer velocity relative to the medium (positive when moving towards source, negative when moving away from source),
## vs is source velocity relative to the medium (positive when moving away from observer, negative when moving towards observer).

# Conditions:
## - Velocities vs and vo are much smaller than v.

observed_frequency = Symbol("observed_frequency", units.frequency)
real_frequency = Symbol("real_frequency", units.frequency)
wave_velocity = Symbol("wave_velocity", units.velocity)
source_velocity = Symbol("source_velocity", units.velocity)
observer_velocity = Symbol("observer_velocity", units.velocity)

law = Eq(observed_frequency,
    real_frequency * (wave_velocity + observer_velocity) / (wave_velocity + source_velocity))


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

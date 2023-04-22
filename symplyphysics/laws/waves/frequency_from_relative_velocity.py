from sympy import (Eq, solve, sqrt)
from sympy.physics.units import speed_of_light
from symplyphysics import (units, expr_to_quantity, Quantity, Symbol,
                           print_expression, validate_input_symbols,
                           validate_output_symbol)

# Description
## Doppler effect is also applicable to electromagnetic waves in vacuum. As there is no any medium required for these waves to propagate,
## speed of source related to observer is used for the calculation of the Doppler effect.

# Law: fo = fs * sqrt((c - v)/(c + v)), where
## fo is observed frequency,
## fs is source wave frequency,
## v is velocity of source related to observer (positive velocity direction means source moving towards observer),
## c is speed of light.

observed_frequency = Symbol("observed_frequency", units.frequency)
real_frequency = Symbol("real_frequency", units.frequency)
source_velocity = Symbol("source_velocity", units.velocity)

law = Eq(
    observed_frequency,
    real_frequency * sqrt((speed_of_light - source_velocity) /
                          (speed_of_light + source_velocity)))


def print() -> str:
    return print_expression(law)


@validate_input_symbols(real_frequency_=real_frequency,
                        source_velocity_=source_velocity)
@validate_output_symbol(observed_frequency)
def calculate_observed_frequency(real_frequency_: Quantity,
                                 source_velocity_: Quantity) -> Quantity:
    result_expr = solve(law, observed_frequency,
                        dict=True)[0][observed_frequency]
    frequency_applied = result_expr.subs({
        real_frequency: real_frequency_,
        source_velocity: source_velocity_
    })
    return expr_to_quantity(frequency_applied)

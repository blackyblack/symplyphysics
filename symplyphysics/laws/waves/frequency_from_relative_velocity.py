from sympy import Expr
from symplyphysics import (
    Eq, pretty, solve, units, expr_to_quantity, sqrt
)
from symplyphysics.core.quantity_decorator import validate_input_symbols, validate_output_symbol
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.core.symbols.symbols import Symbol, to_printable
from sympy.physics.units import speed_of_light

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

law = Eq(observed_frequency, real_frequency * sqrt((speed_of_light - source_velocity) / (speed_of_light + source_velocity)))

def print(expr: Expr) -> str:
    symbols = [observed_frequency, real_frequency, source_velocity]
    return pretty(to_printable(expr, symbols), use_unicode=False)

@validate_input_symbols(real_frequency_=real_frequency, source_velocity_=source_velocity)
@validate_output_symbol(observed_frequency)
def calculate_observed_frequency(real_frequency_: Quantity, source_velocity_: Quantity) -> Quantity:        
    result_expr = solve(law, observed_frequency, dict=True)[0][observed_frequency]
    frequency_applied = result_expr.subs({
        real_frequency: real_frequency_,        
        source_velocity: source_velocity_})
    return expr_to_quantity(frequency_applied)

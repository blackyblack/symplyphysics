from sympy import (Eq, solve, pi)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, dimensionless, convert_to_float)

# Description
## The quality factor shows the ratio of the energy stored in the resonator to the energy loss during one oscillation period.
## If the resonator is an oscillatory circuit, that quality factor will depend on the resistance, inductance, and oscillation frequency.

## Law is: Q = R / (2 * pi * f * L), where
## Q - quality factor of the resonator,
## R - resistance in the oscillating circuit,
## L - inductance in the oscillatory circuit,
## f - oscillation frequency.

quality_factor = Symbol("quality_factor", dimensionless)

resistance = Symbol("resistance", units.impedance)
inductance = Symbol("inductance", units.inductance)
frequency = Symbol("frequency", units.frequency)

law = Eq(quality_factor, resistance / (2 * pi * frequency * inductance))


def print_law() -> str:
    return print_expression(law)


@validate_input(resistance_=resistance, inductance_=inductance, frequency_=frequency)
@validate_output(quality_factor)
def calculate_quality_factor(resistance_: Quantity, inductance_: Quantity,
    frequency_: Quantity) -> float:
    result_expr = solve(law, quality_factor, dict=True)[0][quality_factor]
    result_expr = result_expr.subs({
        resistance: resistance_,
        inductance: inductance_,
        frequency: frequency_,
    })
    return convert_to_float(result_expr)

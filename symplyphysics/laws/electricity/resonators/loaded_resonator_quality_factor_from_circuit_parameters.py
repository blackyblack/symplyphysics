from sympy import (Eq, solve, pi)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input, validate_output,
    dimensionless, convert_to_float)

# Description
## The quality factor shows the ratio of the energy stored in the resonator to the energy loss during one oscillation period.
## If the resonator is an oscillatory circuit to which an external circuit is connected, then the loaded Q-factor of the resonator
## depends on the resistances of the resonator and the external circuit, as well as on the inductance of the resonator and the oscillation frequency.

## Law is: Q = (Rl * R) / (2 * pi * f * L * (Rl + R)), where
## Q - quality factor of the resonator,
## R - resistance in the oscillating circuit,
## L - inductance in the oscillatory circuit,
## f - oscillation frequency,
## Rl - load resistance.

loaded_resonator_quality_factor = Symbol("loaded_resonator_quality_factor", dimensionless)

resistance = Symbol("resistance", units.impedance)
inductance = Symbol("inductance", units.inductance)
frequency = Symbol("frequency", units.frequency)
load_resistance = Symbol("load_resistance", units.impedance)

law = Eq(loaded_resonator_quality_factor, (load_resistance * resistance) / (2 * pi * frequency * inductance * (load_resistance + resistance)))


def print_law() -> str:
    return print_expression(law)


@validate_input(resistance_=resistance,
    inductance_=inductance,
    frequency_=frequency,
    load_resistance_=load_resistance)
@validate_output(loaded_resonator_quality_factor)
def calculate_quality_factor(resistance_: Quantity, inductance_: Quantity, frequency_: Quantity, load_resistance_: Quantity) -> float:
    result_expr = solve(law, loaded_resonator_quality_factor, dict=True)[0][loaded_resonator_quality_factor]
    result_expr = result_expr.subs({
        resistance: resistance_,
        inductance: inductance_,
        frequency: frequency_,
        load_resistance: load_resistance_,
    })
    return convert_to_float(result_expr)

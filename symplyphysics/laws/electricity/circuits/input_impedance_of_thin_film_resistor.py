from sympy import Eq, solve, I, pi
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output)

# Description
## Thin-film resistors in integrated design are used in microwave circuits. The input impedance of a thin-film resistor
## depends on its resistance and capacitance, as well as on the frequency.

## Law is: Zin = R / (1 + I * 2 * pi * f * R * C / 3), where
## Zin - input impedance of a thin-film resistor,
## R - resistance,
## I - imaginary unit,
## f - frequency,
## C - capacitance.

input_impedance = Symbol("input_impedance", units.impedance)

resistance = Symbol("resistance", units.impedance)
frequency = Symbol("frequency", units.frequency)
capacitance = Symbol("capacitance", units.capacitance)

law = Eq(input_impedance, resistance / (1 + I * 2 * pi * frequency * resistance * capacitance / 3))


@validate_input(
    resistance_=resistance,
    frequency_=frequency,
    capacitance_=capacitance,
)
@validate_output(input_impedance)
def calculate_input_impedance(
    resistance_: Quantity,
    frequency_: Quantity,
    capacitance_: Quantity,
) -> Quantity:
    result_expr = solve(law, input_impedance, dict=True)[0][input_impedance]
    result_expr = result_expr.subs({
        resistance: resistance_,
        frequency: frequency_,
        capacitance: capacitance_,
    })
    return Quantity(result_expr)

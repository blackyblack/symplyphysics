from sympy import Eq, solve
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, dimensionless, convert_to_float)

# Description
## Standing wave coefficient is the ratio of the highest voltage in the transmission line to the lowest voltage. It is a measure of matching the load 
## with the transmission line. The coefficient in the transmission line does not depend on the internal resistance of the electromagnetic
## wave source and (in the case of a linear load) on the power of the generator.

## Law is: K = Umax / Umin, where
## K - standing wave coefficient,
## Umax - maximum voltage in the transmission line,
## Umin - minimum voltage in the transmission line.

coefficient_standing_wave = Symbol("coefficient_standing_wave", dimensionless)

maximum_voltage = Symbol("maximum_voltage", units.voltage)
minimum_voltage = Symbol("minimum_voltage", units.voltage)

law = Eq(coefficient_standing_wave, maximum_voltage / minimum_voltage)


def print_law() -> str:
    return print_expression(law)


@validate_input(maximum_voltage_=maximum_voltage,
    minimum_voltage_=minimum_voltage)
@validate_output(coefficient_standing_wave)
def calculate_coefficient_standing_wave(maximum_voltage_: Quantity,
    minimum_voltage_: Quantity) -> float:
    result_expr = solve(law, coefficient_standing_wave, dict=True)[0][coefficient_standing_wave]
    result_expr = result_expr.subs({
        maximum_voltage: maximum_voltage_,
        minimum_voltage: minimum_voltage_,
    })
    return convert_to_float(result_expr)

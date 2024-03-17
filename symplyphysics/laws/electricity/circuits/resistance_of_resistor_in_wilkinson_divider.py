from sympy import (Eq, solve,)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    dimensionless
)

# Description
## The Wilkinson divider is a device designed to divide the power of a microwave signal into two output ports.
## The microstrip version has one surface-mounted resistor. The resistance of this resistor depends on the resistance
## of the transmission line and the power ratio at the output ports.

## Law is: R = Z0 * (1 + k^2) / k, where
## R - resistance of resistor,
## Z0 - transmission line resistance,
## k - ratio coefficient of the power at the outputs of the divider.

resistance = Symbol("resistance", units.impedance)

transmission_line_resistance = Symbol("transmission_line_resistance", units.impedance)
ratio_of_power = Symbol("ratio_of_power", dimensionless)

law = Eq(resistance, transmission_line_resistance * (1 + ratio_of_power**2) / ratio_of_power)


def print_law() -> str:
    return print_expression(law)


@validate_input(transmission_line_resistance_=transmission_line_resistance, ratio_of_power_=ratio_of_power)
@validate_output(resistance)
def calculate_resistance(transmission_line_resistance_: Quantity, ratio_of_power_: float) -> Quantity:
    result_expr = solve(law, resistance, dict=True)[0][resistance]
    result_expr = result_expr.subs({
        transmission_line_resistance: transmission_line_resistance_,
        ratio_of_power: ratio_of_power_,
    })
    return Quantity(result_expr)

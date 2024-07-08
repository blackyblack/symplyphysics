from sympy import Eq, Rational, solve, sqrt
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

# Description
## The steepness of the volt-ampere characteristic is generally equal to the voltage derivative of the current.
## For a vacuum diode, it can be calculated if the diode parameter and the voltage between the anode and cathode are known.

## Law is: S = (3 / 2) * g * sqrt(Ua), where
## S - the steepness of the volt-ampere characteristic,
## g - diode constant,
## Ua - voltage between cathode and anode.

steepness = Symbol("steepness", units.current / units.voltage)

diode_constant = Symbol("diode_constant", units.current / units.voltage**Rational(3, 2))
voltage = Symbol("voltage", units.voltage)

law = Eq(steepness, (3 / 2) * diode_constant * sqrt(voltage))


@validate_input(diode_constant_=diode_constant,
    voltage_=voltage)
@validate_output(steepness)
def calculate_steepness(diode_constant_: Quantity, voltage_: Quantity) -> Quantity:
    result_expr = solve(law, steepness,
        dict=True)[0][steepness]
    result_expr = result_expr.subs({
        diode_constant: diode_constant_,
        voltage: voltage_,
    })
    return Quantity(result_expr)

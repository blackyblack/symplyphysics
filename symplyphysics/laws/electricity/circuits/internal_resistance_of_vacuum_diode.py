from sympy import Eq, solve, sqrt
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

# Description
## The internal resistance is generally equal to the current derivative of the voltage.
## For a vacuum diode, it can be calculated if the diode parameter and the voltage between the anode and cathode are known.

## Law is: ri = 2 / (3 * g * sqrt(Ua)), where
## ri - the internal resistance of the vacuum diode,
## g - diode constant,
## Ua - voltage between cathode and anode.

internal_resistance = Symbol("internal_resistance", units.impedance)

diode_constant = Symbol("diode_constant", units.current / units.voltage**(3 / 2))
voltage = Symbol("voltage", units.voltage)

law = Eq(internal_resistance, 2 / (3 * diode_constant * sqrt(voltage)))


@validate_input(diode_constant_=diode_constant,
    voltage_=voltage)
@validate_output(internal_resistance)
def calculate_internal_resistance(diode_constant_: Quantity, voltage_: Quantity) -> Quantity:
    result_expr = solve(law, internal_resistance,
        dict=True)[0][internal_resistance]
    result_expr = result_expr.subs({
        diode_constant: diode_constant_,
        voltage: voltage_,
    })
    return Quantity(result_expr)

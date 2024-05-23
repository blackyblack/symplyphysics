from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

# Description
## The current-voltage characteristic of a vacuum diode is described by the 3/2-power law. The diode
## constant in this law depends only on the relative position, shape and size of the electrodes of the vacuum diode.

## Law is: I = g * U^(3 / 2), where
## I - anode current,
## g - diode constant,
## U - voltage between cathode and anode.

current = Symbol("current", units.current)

diode_constant = Symbol("diode_constant", units.current / units.voltage**(3 / 2))
voltage = Symbol("voltage", units.voltage)

law = Eq(current, diode_constant * voltage**(3 / 2))


@validate_input(diode_constant_=diode_constant,
    voltage_=voltage)
@validate_output(current)
def calculate_current(diode_constant_: Quantity, voltage_: Quantity) -> Quantity:
    result_expr = solve(law, current,
        dict=True)[0][current]
    result_expr = result_expr.subs({
        diode_constant: diode_constant_,
        voltage: voltage_,
    })
    return Quantity(result_expr)

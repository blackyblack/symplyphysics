from sympy import Eq, Rational, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    dimensionless
)

# Description
## A triode has three electrodes: a cathode, an anode and one control grid. The triode can be replaced with an equivalent diode
## and the 3/2-power law can be applied. The triode constant in this law depends only on the relative position, shape and size of
## the electrodes of the vacuum triode.

## Law is: I = g * (Ua + mu * Uc)^(3 / 2) , where
## I - anode current,
## g - triode constant,
## Ua - voltage between cathode and anode.
## mu - voltage triode gain,
## Uc - grid voltage.

current = Symbol("current", units.current)

triode_constant = Symbol("triode_constant", units.current / units.voltage**Rational(3, 2))
anode_voltage = Symbol("anode_voltage", units.voltage)
voltage_triode_gain = Symbol("voltage_triode_gain", dimensionless)
grid_voltage = Symbol("grid_voltage", units.voltage)

law = Eq(current, triode_constant * (anode_voltage + voltage_triode_gain * grid_voltage)**Rational(3, 2))


@validate_input(triode_constant_=triode_constant,
    anode_voltage_=anode_voltage,
    voltage_triode_gain_=voltage_triode_gain,
    grid_voltage_=grid_voltage)
@validate_output(current)
def calculate_current(triode_constant_: Quantity, anode_voltage_: Quantity,
                      voltage_triode_gain_: float, grid_voltage_: Quantity) -> Quantity:
    result_expr = solve(law, current,
        dict=True)[0][current]
    result_expr = result_expr.subs({
        triode_constant: triode_constant_,
        anode_voltage: anode_voltage_,
        voltage_triode_gain: voltage_triode_gain_,
        grid_voltage: grid_voltage_,
    })
    return Quantity(result_expr)

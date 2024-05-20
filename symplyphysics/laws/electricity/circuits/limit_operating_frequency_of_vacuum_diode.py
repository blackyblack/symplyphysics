from sympy import Eq, solve, sqrt
from sympy.physics.units import elementary_charge, electron_rest_mass
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

# Description
## The vacuum diode has a limit operating frequency. At frequencies higher, the diode stops functioning normally, since the electrons
## will not have time to move from the cathode to the anode.

## Law is: f = sqrt(2 * e * Ua / m) / (6 * d), where
## f - the limit operating frequency of the vacuum diode,
## e - elementary charge,
## Ua - voltage between cathode and anode,
## m - electron rest mass,
## d - distance between electrodes.

limit_operating_frequency = Symbol("limit_operating_frequency", units.frequency)

distance_between_electrodes = Symbol("distance_between_electrodes", units.length)
voltage = Symbol("voltage", units.voltage)

law = Eq(limit_operating_frequency, sqrt(2 * elementary_charge * voltage / electron_rest_mass) / (6 * distance_between_electrodes))


@validate_input(distance_between_electrodes_=distance_between_electrodes,
    voltage_=voltage)
@validate_output(limit_operating_frequency)
def calculate_limit_operating_frequency(distance_between_electrodes_: Quantity, voltage_: Quantity) -> Quantity:
    result_expr = solve(law, limit_operating_frequency,
        dict=True)[0][limit_operating_frequency]
    result_expr = result_expr.subs({
        distance_between_electrodes: distance_between_electrodes_,
        voltage: voltage_,
    })
    return Quantity(result_expr)

from sympy import Eq, solve
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
## and the 3/2-power law can be applied. The voltage in the equivalent diode can also be calculated.

## Law is: U = (Ug + Ua / mu) / (1 + ((Xa / Xg)^(4 / 3)) / mu), where
## U - the voltage between the cathode and the anode of an equivalent diode for a triode,
## Ug - voltage between cathode and grid,
## Ua - voltage between cathode and anode,
## mu - voltage triode gain,
## Xa - the distance from the cathode to the anode,
## Xg - the distance from the cathode to the grid.

voltage_of_equivalent_diode = Symbol("voltage_of_equivalent_diode", units.voltage)

distance_to_anode = Symbol("distance_to_anode", units.length)
distance_to_grid = Symbol("distance_to_grid", units.length)
anode_voltage = Symbol("anode_voltage", units.voltage)
voltage_triode_gain = Symbol("voltage_triode_gain", dimensionless)
grid_voltage = Symbol("grid_voltage", units.voltage)

law = Eq(voltage_of_equivalent_diode, (grid_voltage + anode_voltage / voltage_triode_gain) / (1 + ((distance_to_anode / distance_to_grid)**(4 / 3)) / voltage_triode_gain))


@validate_input(distance_to_anode_=distance_to_anode,
    distance_to_grid_=distance_to_grid,
    anode_voltage_=anode_voltage,
    voltage_triode_gain_=voltage_triode_gain,
    grid_voltage_=grid_voltage)
@validate_output(voltage_of_equivalent_diode)
def calculate_voltage_of_equivalent_diode(distance_to_anode_: Quantity, distance_to_grid_: Quantity, anode_voltage_: Quantity,
                      voltage_triode_gain_: float, grid_voltage_: Quantity) -> Quantity:
    result_expr = solve(law, voltage_of_equivalent_diode,
        dict=True)[0][voltage_of_equivalent_diode]
    result_expr = result_expr.subs({
        distance_to_anode: distance_to_anode_,
        distance_to_grid: distance_to_grid_,
        anode_voltage: anode_voltage_,
        voltage_triode_gain: voltage_triode_gain_,
        grid_voltage: grid_voltage_,
    })
    return Quantity(result_expr)

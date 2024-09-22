from sympy import Eq, Rational, solve
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output, dimensionless)

# Description
## A pentode has five electrodes: a cathode, an anode, and three grids. The pentode can be replaced with an equivalent diode,
## and the voltage in the equivalent diode can be calculated.
## The coefficient of direct permeability of the grid characterizes its shielding effect and shows how much the electrostatic field of the anode is
## weaker than the grid field affects the cathode area.

## Law is: U = (Ug1 + D1 * Ug2 + D1 * D2 * Ug3 + D1 * D2 * D3 * Ua) / (1 + ((Xa / Xg1)^(4 / 3)) * D1), where
## U - the voltage between the cathode and the anode of an equivalent diode for a pentode,
## Ug1 - voltage between cathode and first grid,
## Ug2 - voltage between cathode and second grid,
## Ug3 - voltage between cathode and third grid,
## Ua - voltage between cathode and anode.
## D1 - coefficient of direct permeability of the first grid,
## D2 - coefficient of direct permeability of the second grid,
## D3 - coefficient of direct permeability of the third grid,
## Xa - the distance from the cathode to the anode,
## Xg1 - the distance from the cathode to the first grid.

voltage_of_equivalent_diode = Symbol("voltage_of_equivalent_diode", units.voltage)

voltage_of_first_grid = Symbol("voltage_of_first_grid", units.voltage)
voltage_of_second_grid = Symbol("voltage_of_second_grid", units.voltage)
voltage_of_third_grid = Symbol("voltage_of_third_grid", units.voltage)
anode_voltage = Symbol("anode_voltage", units.voltage)
coefficient_direct_permeability_of_first_grid = Symbol(
    "coefficient_direct_permeability_of_first_grid", dimensionless)
coefficient_direct_permeability_of_second_grid = Symbol(
    "coefficient_direct_permeability_of_second_grid", dimensionless)
coefficient_direct_permeability_of_third_grid = Symbol(
    "coefficient_direct_permeability_of_third_grid", dimensionless)
distance_to_anode = Symbol("distance_to_anode", units.length)
distance_to_first_grid = Symbol("distance_to_first_grid", units.length)

expression_1 = voltage_of_first_grid + coefficient_direct_permeability_of_first_grid * voltage_of_second_grid + coefficient_direct_permeability_of_first_grid * coefficient_direct_permeability_of_second_grid * voltage_of_third_grid
expression_2 = coefficient_direct_permeability_of_first_grid * coefficient_direct_permeability_of_second_grid * coefficient_direct_permeability_of_third_grid * anode_voltage
expression_3 = 1 + ((distance_to_anode / distance_to_first_grid)**Rational(4,
    3)) * coefficient_direct_permeability_of_first_grid
law = Eq(voltage_of_equivalent_diode, (expression_1 + expression_2) / expression_3)


@validate_input(voltage_of_first_grid_=voltage_of_first_grid,
    voltage_of_second_grid_=voltage_of_second_grid,
    voltage_of_third_grid_=voltage_of_third_grid,
    anode_voltage_=anode_voltage,
    coefficient_direct_permeability_of_first_grid_=coefficient_direct_permeability_of_first_grid,
    coefficient_direct_permeability_of_second_grid_=coefficient_direct_permeability_of_second_grid,
    coefficient_direct_permeability_of_third_grid_=coefficient_direct_permeability_of_third_grid,
    distance_to_anode_=distance_to_anode,
    distance_to_first_grid_=distance_to_first_grid)
@validate_output(voltage_of_equivalent_diode)
def calculate_voltage_of_equivalent_diode(voltage_of_first_grid_: Quantity,
    voltage_of_second_grid_: Quantity, voltage_of_third_grid_: Quantity, anode_voltage_: Quantity,
    coefficient_direct_permeability_of_first_grid_: float,
    coefficient_direct_permeability_of_second_grid_: float,
    coefficient_direct_permeability_of_third_grid_: float, distance_to_anode_: Quantity,
    distance_to_first_grid_: Quantity) -> Quantity:
    # pylint: disable=too-many-arguments, too-many-positional-arguments
    if distance_to_anode_.scale_factor <= distance_to_first_grid_.scale_factor:
        raise ValueError(
            "The distance to the anode should be greater than the distance to the grid.")
    result_expr = solve(law, voltage_of_equivalent_diode, dict=True)[0][voltage_of_equivalent_diode]
    result_expr = result_expr.subs({
        voltage_of_first_grid:
            voltage_of_first_grid_,
        voltage_of_second_grid:
            voltage_of_second_grid_,
        voltage_of_third_grid:
            voltage_of_third_grid_,
        anode_voltage:
            anode_voltage_,
        coefficient_direct_permeability_of_first_grid:
            coefficient_direct_permeability_of_first_grid_,
        coefficient_direct_permeability_of_second_grid:
            coefficient_direct_permeability_of_second_grid_,
        coefficient_direct_permeability_of_third_grid:
            coefficient_direct_permeability_of_third_grid_,
        distance_to_anode:
            distance_to_anode_,
        distance_to_first_grid:
            distance_to_first_grid_
    })
    return Quantity(result_expr)

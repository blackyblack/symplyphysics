from sympy import Eq, solve
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output, dimensionless,
    convert_to_float)

# Description
## The coefficient of direct permeability of the grid characterizes its shielding effect and shows how much the electrostatic field of the anode is
## weaker than the grid field affects the cathode area.
## There is no general formula for this coefficient for an arbitrary electrode configuration, but there are formulas for special cases.
## One such case is a triode with flat electrodes.

## Law is: D = T * h / (x * d), where
## D - direct permeability coefficient of the grid,
## T - first tabular coefficient,
## h - grid step(distance between the grid wires),
## x - distance from anode to grid,
## d - second tabular coefficient.

direct_permeability_coefficient = Symbol("direct_permeability_coefficient", dimensionless)

tabular_coefficient_T = Symbol("tabular_coefficient_T", dimensionless)
grid_step = Symbol("grid_step", units.length)
distance_to_grid = Symbol("distance_to_grid", units.length)
tabular_coefficient_d = Symbol("tabular_coefficient_d", dimensionless)

law = Eq(direct_permeability_coefficient,
    tabular_coefficient_T * grid_step / (distance_to_grid * tabular_coefficient_d))


@validate_input(tabular_coefficient_T_=tabular_coefficient_T,
    grid_step_=grid_step,
    distance_to_grid_=distance_to_grid,
    tabular_coefficient_d_=tabular_coefficient_d)
@validate_output(direct_permeability_coefficient)
def calculate_direct_permeability_coefficient(tabular_coefficient_T_: float, grid_step_: Quantity,
    distance_to_grid_: Quantity, tabular_coefficient_d_: float) -> float:
    result_expr = solve(law, direct_permeability_coefficient,
        dict=True)[0][direct_permeability_coefficient]
    result_expr = result_expr.subs({
        tabular_coefficient_T: tabular_coefficient_T_,
        grid_step: grid_step_,
        distance_to_grid: distance_to_grid_,
        tabular_coefficient_d: tabular_coefficient_d_
    })
    return convert_to_float(Quantity(result_expr))

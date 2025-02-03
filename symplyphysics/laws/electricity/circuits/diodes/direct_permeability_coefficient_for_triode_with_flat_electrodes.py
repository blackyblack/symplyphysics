"""
Direct permeability coefficient of triode with flat electrodes
==============================================================

The coefficient of direct permeability of the grid characterizes its shielding effect
and shows how much the electrostatic field of the anode is weaker than the grid field
affects the cathode area. There is no general formula for this coefficient for an
arbitrary electrode configuration, but there are formulas for special cases. One such
case is a triode with flat electrodes.

..
    TODO: find link
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
    dimensionless,
    convert_to_float,
    symbols,
    clone_as_symbol,
)

direct_permeability_coefficient = SymbolNew("D", dimensionless)
"""
Direct permeability coefficient of the grid.
"""

first_tabular_coefficient = SymbolNew("C_1", dimensionless)
"""
First tabular coefficient.
"""

grid_step = clone_as_symbol(symbols.euclidean_distance, subscript="0")
"""
Grid step, or :symbols:`euclidean_distance` between grid wires.
"""

distance = symbols.euclidean_distance
"""
:symbols:`euclidean_distance` between anode and grid.
"""

second_tabular_coefficient = SymbolNew("C_2", dimensionless)
"""
Second tabular coefficient.
"""

law = Eq(direct_permeability_coefficient,
    first_tabular_coefficient * grid_step / (distance * second_tabular_coefficient))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(tabular_coefficient_t_=first_tabular_coefficient,
    grid_step_=grid_step,
    distance_to_grid_=distance,
    tabular_coefficient_d_=second_tabular_coefficient)
@validate_output(direct_permeability_coefficient)
def calculate_direct_permeability_coefficient(tabular_coefficient_t_: float, grid_step_: Quantity,
    distance_to_grid_: Quantity, tabular_coefficient_d_: float) -> float:
    result_expr = solve(law, direct_permeability_coefficient,
        dict=True)[0][direct_permeability_coefficient]
    result_expr = result_expr.subs({
        first_tabular_coefficient: tabular_coefficient_t_,
        grid_step: grid_step_,
        distance: distance_to_grid_,
        second_tabular_coefficient: tabular_coefficient_d_
    })
    return convert_to_float(Quantity(result_expr))

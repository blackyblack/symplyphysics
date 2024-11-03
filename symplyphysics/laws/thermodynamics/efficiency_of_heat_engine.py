"""
Efficiency of heat engine
=========================

Efficiency of a heat engine is the ratio of the useful energy to the total energy received
by the system.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    convert_to_float,
    symbols,
    clone_as_symbol,
)

heat_from_heater = clone_as_symbol(symbols.heat, subscript="h")
r"""
:symbols:`heat` transferred by the heater to the heat engine.
"""

heat_to_refrigerator = clone_as_symbol(symbols.heat, subscript="r")
r"""
:symbols:`heat` transferred by the heat engine to the refrigerator.
"""

efficiency = symbols.thermal_efficiency
"""
:symbols:`thermal_efficiency` of the heat engine.
"""

law = Eq(efficiency, 1 - heat_to_refrigerator / heat_from_heater)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(heat_from_heater_=heat_from_heater, heat_to_refrigerator_=heat_to_refrigerator)
@validate_output(efficiency)
def calculate_efficiency_factor(heat_from_heater_: Quantity,
    heat_to_refrigerator_: Quantity) -> float:
    result_efficiency_factor = solve(law, efficiency, dict=True)[0][efficiency]
    result_expr = result_efficiency_factor.subs({
        heat_from_heater: heat_from_heater_,
        heat_to_refrigerator: heat_to_refrigerator_
    })
    return convert_to_float(result_expr)

"""
Efficiency of heat engine
=========================

Efficiency of a heat engine is the ratio of the useful energy to the total energy received
by the system.
"""

from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    dimensionless,
    convert_to_float,
)

heat_from_heater = Symbol("heat_from_heater", units.energy)
r"""
Heat transferred by the heater to the heat engine.

Symbol:
    :code:`Q_h`

Latex:
    :math:`Q_\text{h}`
"""

heat_to_refrigerator = Symbol("heat_to_refrigerator", units.energy)
r"""
Heat transferred by the heat engine to the refrigerator.

Symbol:
    :code:`Q_r`

Latex:
    :math:`Q_\text{r}`
"""

efficiency = Symbol("efficiency", dimensionless)
r"""
Efficiency of the heat engine.

Symbol:
    :code:`eta`

Latex:
    :math:`\eta`
"""

law = Eq(efficiency, 1 - heat_to_refrigerator / heat_from_heater)
r"""
:code:`eta = 1 - Q_r / Q_h`

Latex:
    .. math::
        \eta = 1 - \frac{Q_\text{r}}{Q_\text{h}}
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

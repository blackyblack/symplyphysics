"""
Heat is heat capacity times temperature change
==============================================

The amount of heat a body receives or releases when its temperature changes is equal
to the heat capacity of the body times the change in the body's temperature.

**Notes:**

#. Usually intensive heat capacity is known, in that case refer to :doc:`laws.quantities` laws.

**Condition:**

#. The heat capacity must be independent of the body's temperature. If this doesn't hold,
   refer to the differential definition of the heat capacity.
"""

from sympy import (Eq, solve)
from symplyphysics import (clone_symbol, symbols, units, Quantity, Symbol, validate_input,
    validate_output)

heat = Symbol("heat", units.energy)
"""
Heat received or released by the body.

Symbol:
    :code:`Q`
"""

heat_capacity = Symbol("heat_capacity", units.energy / units.temperature)
"""
Heat capacity of the body.

Symbol:
    :code:`C`
"""

temperature_change = clone_symbol(symbols.temperature,
    display_symbol="dT",
    display_latex="\\Delta T")
"""
Change in the body's :attr:`~symplyphysics.symbols.temperature`
"""

law = Eq(heat, heat_capacity * temperature_change)
r"""
:code:`Q = C * dT`

Latex:
    .. math::
        Q = C \Delta T
"""


@validate_input(heat_capacity_=heat_capacity, temperature_change_=temperature_change)
@validate_output(heat)
def calculate_amount_energy(heat_capacity_: Quantity, temperature_change_: Quantity) -> Quantity:
    result_amount_energy_expr = solve(law, heat, dict=True)[0][heat]
    result_expr = result_amount_energy_expr.subs({
        heat_capacity: heat_capacity_,
        temperature_change: temperature_change_,
    })
    return Quantity(result_expr)

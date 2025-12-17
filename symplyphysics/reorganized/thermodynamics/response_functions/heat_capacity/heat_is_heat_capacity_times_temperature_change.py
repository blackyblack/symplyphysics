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

**Links:**

#. `Wikipedia, derivable from here <https://en.wikipedia.org/wiki/Heat_capacity#Basic_definition>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    clone_as_symbol,
    symbols,
    Quantity,
    validate_input,
    validate_output,
)

heat = symbols.heat
"""
:symbols:`heat` received or released by the body.
"""

heat_capacity = symbols.heat_capacity
"""
:symbols:`heat_capacity` of the body.
"""

temperature_change = clone_as_symbol(
    symbols.temperature,
    display_symbol="Delta(T)",
    display_latex="\\Delta T",
)
"""
Change in the body's :symbols:`temperature`
"""

law = Eq(heat, heat_capacity * temperature_change)
"""
:laws:symbol::

:laws:latex::
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

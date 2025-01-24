"""
Heat of combustion via mass
===========================

*Heat of combustion* of the heat released during the complete combustion of a body.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Heat_of_combustion>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

heat = symbols.heat
"""
:symbols:`heat` released during combustion.
"""

mass = symbols.mass
"""
:symbols:`mass` of the body subjected to combustion.
"""

specific_heat_of_combustion = clone_as_symbol(symbols.specific_energy, subscript="q")
"""
Heat of combustion per unit mass of the body. See :symbols:`specific_energy`.
"""

law = Eq(heat, specific_heat_of_combustion * mass)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(specific_heat_combustion_=specific_heat_of_combustion, mass_of_matter_=mass)
@validate_output(heat)
def calculate_amount_energy(specific_heat_combustion_: Quantity,
    mass_of_matter_: Quantity) -> Quantity:
    result_amount_energy_expr = solve(law, heat, dict=True)[0][heat]
    result_expr = result_amount_energy_expr.subs({
        specific_heat_of_combustion: specific_heat_combustion_,
        mass: mass_of_matter_
    })
    return Quantity(result_expr)

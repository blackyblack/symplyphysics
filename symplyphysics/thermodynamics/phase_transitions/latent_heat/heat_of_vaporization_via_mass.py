"""
Heat of vaporization via mass
=============================

*Heat of vaporization* is the heat released during the complete vaporization of a body,
in which it is converted from liquid into gaseous state.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Enthalpy_of_vaporization>`__.
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
:symbols:`heat` released during vaporization.
"""

mass = symbols.mass
"""
:symbols:`mass` of the body subjected to vaporization.
"""

specific_heat_of_vaporization = clone_as_symbol(symbols.specific_energy, subscript="L")
"""
Heat of vaporization per unit mass of the body. See :symbols:`specific_energy`.
"""

law = Eq(heat, specific_heat_of_vaporization * mass)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(specific_heat_vaporization_=specific_heat_of_vaporization, mass_of_matter_=mass)
@validate_output(heat)
def calculate_amount_energy(specific_heat_vaporization_: Quantity,
    mass_of_matter_: Quantity) -> Quantity:
    result_amount_energy_expr = solve(law, heat, dict=True)[0][heat]
    result_expr = result_amount_energy_expr.subs({
        specific_heat_of_vaporization: specific_heat_vaporization_,
        mass: mass_of_matter_
    })
    return Quantity(result_expr)

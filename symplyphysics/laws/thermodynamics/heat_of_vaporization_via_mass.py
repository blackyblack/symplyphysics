"""
Heat of vaporization via mass
=============================

*Heat of vaporization* is the heat released during the complete vaporization of a body,
in which it is converted from liquid into gaseous state.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output, symbols)

heat = Symbol("heat", units.energy)
"""
Heat released during vaporization.

Symbol:
    :code:`Q`
"""

mass = symbols.mass
"""
:attr:`~symplyphysics.symbols.mass` of the body subjected to vaporization.
"""

specific_heat_of_vaporization = Symbol("specific_heat_of_vaporization", units.energy / units.mass)
"""
Heat of vaporization per unit mass of the body.

Symbol:
    :code:`L`
"""

law = Eq(heat, specific_heat_of_vaporization * mass)
r"""
:code:`Q = L * m`

Latex:
    .. math::
        Q = L m
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

"""
Heat of combustion via mass
===========================

*Heat of combustion* of the heat released during the complete combustion of a body.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output, symbols)

heat = Symbol("heat", units.energy)
"""
Heat released during combustion.

Symbol:
    :code:`Q`
"""

mass = symbols.mass
"""
:attr:`~symplyphysics.symbols.mass` of the body subjected to combustion.
"""

specific_heat_of_combustion = Symbol("specific_heat_of_combustion", units.energy / units.mass)
"""
Heat of combustion per unit mass of the body.

Symbol:
    :code:`q`
"""

law = Eq(heat, specific_heat_of_combustion * mass)
"""
:code:`Q = q * m`

Latex:
    .. math::
        Q = q m
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

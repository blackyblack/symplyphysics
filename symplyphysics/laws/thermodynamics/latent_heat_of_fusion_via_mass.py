"""
Latent heat of fusion via mass
==============================

*Latent heat of fusion* is the heat released into or withdrawn from the enviroment when
the body changes its state from a solid to a liquid. The same law applies to the process
of solidification, when the body changes its state from a liquid to a solid, and the
specific heat of solidification is equal by magnitude to that of fusion.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output, symbols)

latent_heat = Symbol("latent_heat", units.energy)
"""
Latent heat of fusion or solidification.

Symbol:
    :code:`Q`
"""

mass = symbols.basic.mass
"""
:attr:`~symplyphysics.symbols.basic.mass` of a melting or solidifying body.

Symbol:
    :code:`m`
"""

specific_heat_of_fusion = Symbol("specific_heat_of_fusion", units.energy / units.mass)
r"""
Heat of fusion or solidification per unit mass.

Symbol:
    :code:`lambda`

Latex:
    :math:`\lambda`
"""

law = Eq(latent_heat, specific_heat_of_fusion * mass)
r"""
:code:`Q = lambda * m`

Latex:
    .. math::
        Q = \lambda m
"""


@validate_input(specific_heat_melting_=specific_heat_of_fusion, mass_of_matter_=mass)
@validate_output(latent_heat)
def calculate_amount_energy(specific_heat_melting_: Quantity,
    mass_of_matter_: Quantity) -> Quantity:
    result_amount_energy_expr = solve(law, latent_heat, dict=True)[0][latent_heat]
    result_expr = result_amount_energy_expr.subs({
        specific_heat_of_fusion: specific_heat_melting_,
        mass: mass_of_matter_
    })
    return Quantity(result_expr)

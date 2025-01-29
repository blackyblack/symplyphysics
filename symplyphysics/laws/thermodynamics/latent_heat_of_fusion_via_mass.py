"""
Latent heat of fusion via mass
==============================

*Latent heat of fusion* is the heat released into or withdrawn from the enviroment when
the body changes its state from a solid to a liquid. The same law applies to the process
of solidification, when the body changes its state from a liquid to a solid, and the
specific heat of solidification is equal by magnitude to that of fusion.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Enthalpy_of_fusion>`__.
"""

from sympy import (Eq, solve)
from symplyphysics import (Quantity, validate_input, validate_output, symbols, clone_as_symbol)

latent_heat = symbols.heat
"""
Latent :symbols:`heat` of fusion or solidification.
"""

mass = symbols.mass
"""
:symbols:`mass` of a melting or solidifying body.
"""

specific_heat_of_fusion = clone_as_symbol(symbols.specific_energy, display_symbol="epsilon_lambda", display_latex="\\varepsilon_\\lambda")
"""
Heat of fusion or solidification per unit mass. See :symbols:`specific_energy`.
"""

law = Eq(latent_heat, specific_heat_of_fusion * mass)
"""
:laws:symbol::

:laws:latex::
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

"""
Diffusion coefficient from macroscopic scattering cross section
===============================================================

Neutron diffusion coefficient can be found from the macroscopic scattering cross
section.

**Links:**

#. `NuclearPower <https://www.nuclear-power.com/nuclear-power/reactor-physics/neutron-diffusion-theory/diffusion-coefficient/>`__.
"""

from sympy import (Eq, solve)
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol
)

macroscopic_transport_cross_section = clone_as_symbol(symbols.macroscopic_cross_section, display_symbol="Sigma_tr", display_latex="\\Sigma_\\text{tr}")
"""
:symbols:`macroscopic_cross_section` of transport.
"""

diffusion_coefficient = symbols.neutron_diffusion_coefficient
"""
:symbols:`neutron_diffusion_coefficient`.
"""

law = Eq(diffusion_coefficient, 1 / (3 * macroscopic_transport_cross_section))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(macroscopic_transport_cross_section_=macroscopic_transport_cross_section)
@validate_output(diffusion_coefficient)
def calculate_diffusion_coefficient(macroscopic_transport_cross_section_: Quantity) -> Quantity:
    result_coefficient_expr = solve(law, diffusion_coefficient,
        dict=True)[0][diffusion_coefficient]
    result_expr = result_coefficient_expr.subs(
        {macroscopic_transport_cross_section: macroscopic_transport_cross_section_})
    return Quantity(result_expr)

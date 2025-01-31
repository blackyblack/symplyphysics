"""
Macroscopic transport cross section from macroscopic scattering cross section
=============================================================================

Macroscopic transport cross section can be found using the macroscopic scattering cross
section and the averaged cosine of the scattering angle.

**Links:**

#. `NuclearPower <https://www.nuclear-power.com/nuclear-power/reactor-physics/neutron-diffusion-theory/diffusion-coefficient/>`__.
"""

from sympy import (Eq, solve)
from symplyphysics import (
    Quantity,
    SymbolNew,
    dimensionless,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

macroscopic_scattering_cross_section = clone_as_symbol(symbols.macroscopic_cross_section, display_symbol="Sigma_s", display_latex="\\sigma_\\text{s}")
"""
:symbols:`macroscopic_cross_section` of scattering.
"""

average_scattering_angle_cosine = SymbolNew("mu", dimensionless, display_latex="\\mu")
"""
Average of the cosine of the angle at which neutrons are scattered in the medium in the
lab system.
"""

macroscopic_transport_cross_section = clone_as_symbol(symbols.macroscopic_cross_section, display_symbol="Sigma_tr", display_latex="\\Sigma_\\text{tr}")
"""
:symbols:`macroscopic_cross_section` of transport.
"""

law = Eq(macroscopic_transport_cross_section,
    macroscopic_scattering_cross_section * (1 - average_scattering_angle_cosine))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(macroscopic_scattering_cross_section_=macroscopic_scattering_cross_section,
    average_scattering_angle_cosine_=average_scattering_angle_cosine)
@validate_output(macroscopic_transport_cross_section)
def calculate_cross_section(macroscopic_scattering_cross_section_: Quantity,
    average_scattering_angle_cosine_: float) -> Quantity:
    result_cross_section_expr = solve(law, macroscopic_transport_cross_section,
        dict=True)[0][macroscopic_transport_cross_section]
    result_expr = result_cross_section_expr.subs({
        macroscopic_scattering_cross_section: macroscopic_scattering_cross_section_,
        average_scattering_angle_cosine: average_scattering_angle_cosine_
    })
    return Quantity(result_expr)

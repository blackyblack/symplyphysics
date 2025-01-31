"""
Diffusion area from diffusion coefficient and absorption cross section
======================================================================

Diffusion area of a nuclear reaction can be found as the ratio of he diffusion
coefficient to the macroscopic absorption cross section of the reaction.

**Links:**

#. `NuclearPower <https://www.nuclear-power.com/nuclear-power/reactor-physics/neutron-diffusion-theory/diffusion-length/>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

diffusion_coefficient = symbols.neutron_diffusion_coefficient
"""
:symbols:`neutron_diffusion_coefficient`.
"""

macroscopic_absorption_cross_section = clone_as_symbol(symbols.macroscopic_cross_section, display_symbol="Sigma_a", display_latex="\\Sigma_\\text{a}")
"""
:symbols:`macroscopic_cross_section` of absorption.
"""

diffusion_area = symbols.neutron_diffusion_area
"""
:symbols:`neutron_diffusion_area`.
"""

law = Eq(diffusion_area, diffusion_coefficient / macroscopic_absorption_cross_section)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(diffusion_coefficient_=diffusion_coefficient,
    macroscopic_absorption_cross_section_=macroscopic_absorption_cross_section)
@validate_output(diffusion_area)
def calculate_diffusion_area(diffusion_coefficient_: Quantity,
    macroscopic_absorption_cross_section_: Quantity) -> Quantity:
    result_diffusion_expr = solve(law, diffusion_area, dict=True)[0][diffusion_area]
    result_expr = result_diffusion_expr.subs({
        diffusion_coefficient: diffusion_coefficient_,
        macroscopic_absorption_cross_section: macroscopic_absorption_cross_section_
    })
    return Quantity(result_expr)

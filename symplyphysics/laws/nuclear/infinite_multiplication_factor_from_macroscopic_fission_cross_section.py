"""
Infinite multiplication factor from macroscopic cross sections
==============================================================

Infinite multiplication factor can be found using the macroscopic fission and absorption
cross section and the average number of neutrons produced per fission.

..
    NOTE: possible link here, but it is about thermal fission factor <https://en.wikipedia.org/wiki/Six_factor_formula#>
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    convert_to_float,
    symbols,
    clone_as_symbol,
)

neutrons_per_fission = clone_as_symbol(symbols.particle_count, display_symbol="nu", display_latex="\\nu")
"""
Average number of neutrons produced per fission. See :symbols:`particle_count`.
"""

macroscopic_fission_cross_section = clone_as_symbol(symbols.macroscopic_cross_section, display_symbol="Sigma_f", display_latex="\\Sigma_text{f}")
"""
:symbols:`macroscopic_cross_section` of fission.
"""

macroscopic_absorption_cross_section = clone_as_symbol(symbols.macroscopic_cross_section, display_symbol="Sigma_a", display_latex="\\Sigma_\\text{a}")
"""
:symbols:`macroscopic_cross_section` of absorption.
"""

infinite_multiplication_factor = symbols.infinite_multiplication_factor
"""
:symbols:`infinite_multiplication_factor`.
"""

law = Eq(
    infinite_multiplication_factor,
    neutrons_per_fission * macroscopic_fission_cross_section / macroscopic_absorption_cross_section)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(neutrons_per_fission_=neutrons_per_fission,
    macroscopic_fission_cross_section_=macroscopic_fission_cross_section,
    macroscopic_absorption_cross_section_=macroscopic_absorption_cross_section)
@validate_output(infinite_multiplication_factor)
def calculate_multiplication_factor(neutrons_per_fission_: float,
    macroscopic_fission_cross_section_: Quantity,
    macroscopic_absorption_cross_section_: Quantity) -> float:

    result_factor_expr = solve(law, infinite_multiplication_factor,
        dict=True)[0][infinite_multiplication_factor]
    result_expr = result_factor_expr.subs({
        neutrons_per_fission: neutrons_per_fission_,
        macroscopic_fission_cross_section: macroscopic_fission_cross_section_,
        macroscopic_absorption_cross_section: macroscopic_absorption_cross_section_
    })
    return convert_to_float(result_expr)

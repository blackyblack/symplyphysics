"""
Reproduction factor from macroscopic cross sections in fuel
===========================================================

Reproduction factor can be found using the macroscopic fissure and absorption cross
sections in the fuel and the average number of neutrons produced per fission.

**Links:**

#. `Wikipedia, first row in table <https://en.wikipedia.org/wiki/Four_factor_formula>`__.
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

macroscopic_fuel_fission_cross_section = clone_as_symbol(
    symbols.macroscopic_cross_section,
    display_symbol="Sigma_ff",
    display_latex="\\Sigma_\\text{f}^\\text{f}",
)
"""
:symbols:`macroscopic_cross_section` of fission in the fuel.
"""

macroscopic_fuel_absorption_cross_section = clone_as_symbol(
    symbols.macroscopic_cross_section,
    display_symbol="Sigma_af",
    display_latex="\\Sigma_\\text{a}^\\text{f}",
)
"""
:symbols:`macroscopic_cross_section` of absorption in the fuel.
"""

reproduction_factor = symbols.reproduction_factor
"""
Neutron :symbols:`reproduction_factor`.
"""

law = Eq(
    reproduction_factor, neutrons_per_fission * macroscopic_fuel_fission_cross_section /
    macroscopic_fuel_absorption_cross_section)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(neutrons_per_fission_=neutrons_per_fission,
    macroscopic_fuel_fission_cross_section_=macroscopic_fuel_fission_cross_section,
    macroscopic_fuel_absorption_cross_section_=macroscopic_fuel_absorption_cross_section)
@validate_output(reproduction_factor)
def calculate_reproduction_factor(neutrons_per_fission_: float,
    macroscopic_fuel_fission_cross_section_: Quantity,
    macroscopic_fuel_absorption_cross_section_: Quantity) -> float:

    if macroscopic_fuel_fission_cross_section_.scale_factor > macroscopic_fuel_absorption_cross_section_.scale_factor:
        raise ValueError(
            f"macroscopic_fuel_fission_cross_section_ ({macroscopic_fuel_fission_cross_section_.scale_factor}) should be <= "
            f"macroscopic_fuel_absorption_cross_section_ ({macroscopic_fuel_absorption_cross_section_.scale_factor})"
        )

    result_factor_expr = solve(law, reproduction_factor,
        dict=True)[0][reproduction_factor]
    result_expr = result_factor_expr.subs({
        neutrons_per_fission: neutrons_per_fission_,
        macroscopic_fuel_fission_cross_section: macroscopic_fuel_fission_cross_section_,
        macroscopic_fuel_absorption_cross_section: macroscopic_fuel_absorption_cross_section_
    })
    return convert_to_float(result_expr)

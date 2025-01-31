"""
Thermal utilization factor from macroscopic absorption cross sections
=====================================================================

Thermal utilization factor can be found using the macroscopic absorption cross section
in the fuel and the total macroscopic absorption cross section in the reactor.

**Links:**

#. `ScienceDirect <https://www.sciencedirect.com/topics/engineering/thermal-utilisation-factor>`__.
#. `NuclearPower <https://www.nuclear-power.com/nuclear-power/reactor-physics/nuclear-fission-chain-reaction/thermal-utilization-factor/>`__.
#. `Wikipedia, second row in table <https://en.wikipedia.org/wiki/Six_factor_formula>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    convert_to_float,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)
from symplyphysics.core.symbols.probability import Probability

macroscopic_fuel_absorption_cross_section = clone_as_symbol(
    symbols.macroscopic_cross_section,
    display_symbol="Sigma_af",
    display_latex="\\Sigma_\\text{a}^\\text{f}",
)
"""
:symbols:`macroscopic_cross_section` of absorption in the fuel.
"""

macroscopic_total_absorption_cross_section = clone_as_symbol(
    symbols.macroscopic_cross_section,
    display_symbol="Sigma_a",
    display_latex="\\Sigma_\\text{a}",
)
"""
Total :symbols:`macroscopic_cross_section` of absorption.
"""

thermal_utilisation_factor = symbols.thermal_utilization_factor
"""
:symbols:`thermal_utilization_factor`.
"""

law = Eq(thermal_utilisation_factor,
    macroscopic_fuel_absorption_cross_section / macroscopic_total_absorption_cross_section)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    macroscopic_fuel_absorption_cross_section_=macroscopic_fuel_absorption_cross_section,
    macroscopic_total_absorption_cross_section_=macroscopic_total_absorption_cross_section)
@validate_output(thermal_utilisation_factor)
def calculate_utilisation_factor(
        macroscopic_fuel_absorption_cross_section_: Quantity,
        macroscopic_total_absorption_cross_section_: Quantity) -> Probability:
    if macroscopic_fuel_absorption_cross_section_.scale_factor > macroscopic_total_absorption_cross_section_.scale_factor:
        raise ValueError(
            f"macroscopic_fuel_absorption_cross_section_ ({macroscopic_fuel_absorption_cross_section_.scale_factor}) should be <= "
            f"macroscopic_total_absorption_cross_section_ ({macroscopic_total_absorption_cross_section_.scale_factor})"
        )

    result_factor_expr = solve(law, thermal_utilisation_factor,
        dict=True)[0][thermal_utilisation_factor]
    result_expr = result_factor_expr.subs({
        macroscopic_fuel_absorption_cross_section: macroscopic_fuel_absorption_cross_section_,
        macroscopic_total_absorption_cross_section: macroscopic_total_absorption_cross_section_
    })
    return Probability(convert_to_float(result_expr))

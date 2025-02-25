"""
Cross section of ionization of atom by electrons per Granovsky
==============================================================

The **Granovsky approximation** for the ionization cross section is considered. As the
electron energy increases, the velocities of primary and secondary electrons increase,
the possibility of their recombination with ions decreases, and the ionization cross
section area increases. However, at very high energies, the ionization cross section
decreases, as the electrons fly past the atom without having time to ionize it, since
the time spent by the electron near the atom decreases.

..
    TODO: find link
    TODO: move to `ionization` folder?
"""

from sympy import Eq, solve, exp
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

effective_cross_section = clone_as_symbol(symbols.cross_section, display_symbol="sigma_eff", display_latex="\\sigma_\\text{eff}")
"""
Effective :symbols:`cross_section` of ionization of particles.
"""

ionization_energy = clone_as_symbol(symbols.energy, display_symbol="sigma_i", display_latex="\\sigma_\\text{i}")
"""
Ionization :symbols:`energy`.
"""

electron_energy = symbols.energy
"""
:symbols:`energy` of ionizing electrons.
"""

maximum_cross_section = clone_as_symbol(symbols.cross_section, display_symbol="sigma_max", display_latex="\\sigma_\\text{max}")
"""
Maximum :symbols:`cross_section` of ionization.
"""

maximum_electron_energy = clone_as_symbol(symbols.energy, display_symbol="E_max", display_latex="E_\\text{max}")
"""
Maximum :symbols:`energy` of ionizing electrons.
"""

law = Eq(
    effective_cross_section,
    maximum_cross_section * ((electron_energy - ionization_energy) /
    (maximum_electron_energy - ionization_energy)) * exp(
    (maximum_electron_energy - electron_energy) / (maximum_electron_energy - ionization_energy)),
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(ionization_energy_=ionization_energy,
    energy_of_electron_=electron_energy,
    maximum_cross_sectional_area_of_ionization_=maximum_cross_section,
    energy_of_electron_at_max_area_=maximum_electron_energy)
@validate_output(effective_cross_section)
def calculate_cross_sectional_area_of_ionization(
        ionization_energy_: Quantity, energy_of_electron_: Quantity,
        maximum_cross_sectional_area_of_ionization_: Quantity,
        energy_of_electron_at_max_area_: Quantity) -> Quantity:
    result_expr = solve(law, effective_cross_section, dict=True)[0][effective_cross_section]
    result_expr = result_expr.subs({
        ionization_energy: ionization_energy_,
        electron_energy: energy_of_electron_,
        maximum_cross_section: maximum_cross_sectional_area_of_ionization_,
        maximum_electron_energy: energy_of_electron_at_max_area_,
    })
    return Quantity(result_expr)

"""
Cross section of ionization of atom by electrons per Granovsky
==============================================================

The effective cross section is a physical quantity characterizing the probability of
transition of a system of two interacting particles to a certain final state, a
quantitative characteristic of the acts of collision of particles of a stream hitting a
target with target particles. The effective cross-section has the dimension of the area.
In this law, we are talking about the interaction of an atom and an electron, which
ionizes an atom. In this case, the Granovsky approximation for the ionization cross
section is considered. As the electron energy increases, the velocities of primary and
secondary electrons increase, the possibility of their recombination with ions decreases,
and the ionization cross-section area increases. However, at very high energies, the
ionization cross section decreases, as the electrons fly past the atom without having
time to ionize it, since the time spent by the electron near the atom decreases.

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

cross_section = symbols.cross_section
"""
:symbols:`cross_section` of ionization of particles.
"""

ionization_energy = clone_as_symbol(symbols.energy, subscript="\\text{i}")
"""
Ionization :symbols:`energy`.
"""

electron_energy = symbols.energy
"""
:symbols:`energy` of ionizing electrons.
"""

maximum_cross_section = clone_as_symbol(symbols.cross_section, subscript="\\text{max}")
"""
Maximum :symbols:`cross_section` of ionization.
"""

maximum_electron_energy = clone_as_symbol(symbols.energy, subscript="\\text{max}")
"""
Maximum :symbols:`energy` of ionizing electrons.
"""

law = Eq(
    cross_section,
    maximum_cross_section
    * ((electron_energy - ionization_energy) / (maximum_electron_energy - ionization_energy))
    * exp((maximum_electron_energy - electron_energy) / (maximum_electron_energy - ionization_energy)),
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(ionization_energy_=ionization_energy,
    energy_of_electron_=electron_energy,
    maximum_cross_sectional_area_of_ionization_=maximum_cross_section,
    energy_of_electron_at_max_area_=maximum_electron_energy)
@validate_output(cross_section)
def calculate_cross_sectional_area_of_ionization(
        ionization_energy_: Quantity, energy_of_electron_: Quantity,
        maximum_cross_sectional_area_of_ionization_: Quantity,
        energy_of_electron_at_max_area_: Quantity) -> Quantity:
    result_expr = solve(law, cross_section,
        dict=True)[0][cross_section]
    result_expr = result_expr.subs({
        ionization_energy: ionization_energy_,
        electron_energy: energy_of_electron_,
        maximum_cross_section: maximum_cross_sectional_area_of_ionization_,
        maximum_electron_energy: energy_of_electron_at_max_area_,
    })
    return Quantity(result_expr)

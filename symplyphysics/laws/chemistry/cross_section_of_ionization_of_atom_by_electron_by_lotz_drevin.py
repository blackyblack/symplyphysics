"""
Cross section of ionization of atom by electrons per Lotz-Drevin
================================================================

In this law, we are talking about the interaction of an atom and an electron, which
ionizes an atom. Equivalent electrons on the outer shell of the ionized atom are
electrons with the same principal and orbital quantum numbers. In this case, the
Lotz-Drewin approximation for the ionization cross section is considered.

**Notation:**

#. :quantity_notation:`bohr_radius`.
#. :quantity_notation:`hydrogen_ionization_energy`.

..
    TODO: find link
    TODO: fix file name
    TODO: move to `ionization` folder?
"""

from sympy import Eq, solve, pi, log
from symplyphysics import (Quantity, Symbol, validate_input, validate_output, dimensionless,
    symbols, clone_as_symbol)
from symplyphysics.quantities import bohr_radius, hydrogen_ionization_energy

cross_section = symbols.cross_section
"""
:symbols:`cross_section` of ionization.
"""

ionization_energy = clone_as_symbol(symbols.energy, subscript="\\text{i}")
"""
:symbols:`energy` of ionization of atoms.
"""

electron_energy = symbols.energy
"""
:symbols:`energy` of ionizing electrons.
"""

first_coefficient = Symbol("A", dimensionless)
"""
A coefficient used in the calculation.
"""

second_coefficient = Symbol("B", dimensionless)
"""
A coefficient used in the calculation.
"""

electron_count = symbols.nonnegative_number
"""
A :symbols:`nonnegative_number` of equivalent electrons on the outer shell of the ionized atom.
"""

_reduced_energy = electron_energy / ionization_energy

law = Eq(
    cross_section,
    (2.66 * pi * bohr_radius**2 * electron_count * hydrogen_ionization_energy**2 /
    ionization_energy**2) * (first_coefficient * (_reduced_energy - 1) / _reduced_energy**2) *
    log(1.25 * second_coefficient * _reduced_energy),
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(ionization_energy_=ionization_energy,
    energy_of_electron_=electron_energy,
    first_calculation_coefficient_=first_coefficient,
    second_calculation_coefficient_=second_coefficient,
    number_of_equivalent_electrons_on_outer_orbit_=electron_count)
@validate_output(cross_section)
def calculate_cross_sectional_area_of_ionization(
        ionization_energy_: Quantity, energy_of_electron_: Quantity,
        first_calculation_coefficient_: float, second_calculation_coefficient_: float,
        number_of_equivalent_electrons_on_outer_orbit_: int) -> Quantity:
    result_expr = solve(law, cross_section, dict=True)[0][cross_section]
    result_expr = result_expr.subs({
        ionization_energy: ionization_energy_,
        electron_energy: energy_of_electron_,
        first_coefficient: first_calculation_coefficient_,
        second_coefficient: second_calculation_coefficient_,
        electron_count: number_of_equivalent_electrons_on_outer_orbit_,
    })
    return Quantity(result_expr)

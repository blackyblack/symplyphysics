"""
Cross section of interaction in Coulomb's interaction model
===========================================================

The effective cross section is a physical quantity characterizing the probability of
transition of a system of two interacting particles to a certain final state, a
quantitative characteristic of the acts of collision of particles of a stream hitting a
target with target particles. The effective cross-section has the dimension of the area.
In a magnetron, this value can be calculated if you know the ionization energy of gas
atoms.

**Notation:**

#. :quantity_notation:`elementary_charge`.
#. :quantity_notation:`vacuum_permittivity`.

..
    TODO: find link
    TODO: move to `magnetron` folder?
"""

from sympy import Eq, solve, pi
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    quantities,
    symbols,
    clone_as_symbol,
)

cross_sectional_area_of_interaction = symbols.cross_section
"""
:symbols:`cross_section` of the interaction of particles.
"""

ionization_energy = clone_as_symbol(symbols.voltage, display_symbol="E_i", display_latex="E_\\text{i}")
"""
Ionization :symbols:`energy` of atoms in terms of :symbols:`voltage`.
"""

law = Eq(cross_sectional_area_of_interaction,
    quantities.elementary_charge**2 / (2 * pi * quantities.vacuum_permittivity**2 * ionization_energy**2))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(ionization_energy_=ionization_energy)
@validate_output(cross_sectional_area_of_interaction)
def calculate_cross_sectional_area_of_interaction(ionization_energy_: Quantity) -> Quantity:
    result_expr = solve(law, cross_sectional_area_of_interaction,
        dict=True)[0][cross_sectional_area_of_interaction]
    result_expr = result_expr.subs({
        ionization_energy: ionization_energy_,
    })
    return Quantity(result_expr)

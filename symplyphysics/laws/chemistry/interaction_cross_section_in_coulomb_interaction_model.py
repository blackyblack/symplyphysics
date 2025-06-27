"""
Interaction cross section in Coulomb's interaction model
========================================================

In a magnetron, the effective cross section of particle interaction can be calculated via the
ionization energy of gas atoms. See :ref:`Effective cross section <_effective_cross_section>`.

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

interaction_cross_section = symbols.cross_section
"""
:symbols:`cross_section` of the interaction of particles.
"""

ionization_energy = clone_as_symbol(symbols.voltage,
    display_symbol="E_i",
    display_latex="E_\\text{i}")
"""
Ionization :symbols:`energy` of atoms expressed in :symbols:`voltage`.
"""

law = Eq(
    interaction_cross_section, quantities.elementary_charge**2 /
    (2 * pi * quantities.vacuum_permittivity**2 * ionization_energy**2))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(ionization_energy_=ionization_energy)
@validate_output(interaction_cross_section)
def calculate_cross_sectional_area_of_interaction(ionization_energy_: Quantity) -> Quantity:
    result_expr = solve(law, interaction_cross_section, dict=True)[0][interaction_cross_section]
    result_expr = result_expr.subs({
        ionization_energy: ionization_energy_,
    })
    return Quantity(result_expr)

"""
Mechanical energy is kinetic and potential energy
=================================================

*Mechanical energy* of the system is defined as the sum of its kinetic energy and potential energy.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

mechanical_energy = symbols.mechanical_energy
"""
:symbols:`mechanical_energy` of the system.
"""

kinetic_energy = symbols.kinetic_energy
"""
:symbols:`kinetic_energy` of the system.
"""

potential_energy = symbols.potential_energy
"""
:symbols:`potential_energy` of the system.
"""

definition = Eq(mechanical_energy, kinetic_energy + potential_energy)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(kinetic_energy_=kinetic_energy, potential_energy_=potential_energy)
@validate_output(mechanical_energy)
def calculate_mechanical_energy(kinetic_energy_: Quantity, potential_energy_: Quantity) -> Quantity:
    solved = solve(definition, mechanical_energy, dict=True)[0][mechanical_energy]
    result_expr = solved.subs({
        kinetic_energy: kinetic_energy_,
        potential_energy: potential_energy_
    })
    return Quantity(result_expr)

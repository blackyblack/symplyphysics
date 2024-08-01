"""
Mechanical energy is kinetic and potential energy
=================================================

*Mechanical energy* of the system is defined as the sum of its kinetic energy and potential energy.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output)

mechanical_energy = Symbol("mechanical_energy", units.energy)
"""
Mechanical energy of the system.

Symbol:
    :code:`E`
"""

kinetic_energy = Symbol("kinetic_energy", units.energy)
"""
Kinetic energy of the system.

Symbol:
    :code:`K`
"""

potential_energy = Symbol("potential_energy", units.energy)
"""
Potential energy of the system.

Symbol:
    :code:`U`
"""

definition = Eq(mechanical_energy, kinetic_energy + potential_energy)
"""
:code:`E = K + U`

Latex:
    .. math::
        E = K + U
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

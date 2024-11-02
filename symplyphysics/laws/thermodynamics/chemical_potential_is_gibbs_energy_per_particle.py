"""
Chemical potential is Gibbs energy per particle
===============================================

Chemical potential of a thermodynamic system can be expressed as the Gibbs energy of the system
per particle. Therefore, chemical potential is an intensive physical quantity, whereas Gibbs energy
and particle count are extensive.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

chemical_potential = symbols.chemical_potential
r"""
:symbols:`chemical_potential` of the system.
"""

gibbs_energy = symbols.gibbs_energy
"""
:symbols:`gibbs_energy` of the system.
"""

particle_count = symbols.particle_count
"""
:symbols:`particle_count` of the system.
"""

law = Eq(chemical_potential, gibbs_energy / particle_count)
"""
:laws:symbol::

:laws:latex::
"""

# TODO: derive law from the extensive property of chemical potential


@validate_input(
    gibbs_energy_=gibbs_energy,
    particle_count_=particle_count,
)
@validate_output(chemical_potential)
def calculate_chemical_potential(
    gibbs_energy_: Quantity,
    particle_count_: int,
) -> Quantity:
    result = law.rhs.subs({
        gibbs_energy: gibbs_energy_,
        particle_count: particle_count_,
    })
    return Quantity(result)

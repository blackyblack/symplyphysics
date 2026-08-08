"""
Chemical potential is Gibbs energy per particle
===============================================

Chemical potential of a thermodynamic system can be expressed as the Gibbs energy of the system
per particle. Therefore, chemical potential is an intensive physical quantity, whereas Gibbs energy
and particle count are extensive.

**Links:**

#. `Wikipedia, last formula in paragraph <https://en.wikipedia.org/wiki/Gibbs_free_energy#Definitions>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.thermodynamics.thermodynamic_potentials.euler_relations import gibbs_energy_euler_relation

chemical_potential = symbols.chemical_potential
"""
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

# Derive the law from the Euler relation for the Gibbs energy.

## The Euler relation states that the Gibbs energy is the product of the chemical potential
## and the particle count, so the chemical potential is found by solving it.
_chemical_potential_derived = solve(
    gibbs_energy_euler_relation.law.subs({
    gibbs_energy_euler_relation.gibbs_energy: gibbs_energy,
    gibbs_energy_euler_relation.chemical_potential: chemical_potential,
    gibbs_energy_euler_relation.particle_count: particle_count,
    }),
    chemical_potential,
)[0]

assert expr_equals(_chemical_potential_derived, law.rhs)


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

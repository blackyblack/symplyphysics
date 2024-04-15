from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## Chemical potential of a thermodynamic system can be viewed as the Gibbs energy of the system
## per particle. Therefore, chemical potential is an intensive physical quantity, whereas Gibbs energy
## and particle count are extensive.

# Law: mu = G / N
## mu - chemical potential
## G - [Gibbs energy](./isobaric_reaction_potential.py)
## N - particle count

chemical_potential = Symbol("chemical_potential", units.energy)
gibbs_energy = Symbol("gibbs_energy", units.energy)
particle_count = Symbol("particle_count")

law = Eq(chemical_potential, gibbs_energy / particle_count)

# TODO: derive law from the extensive property of chemical potential


def print_law() -> str:
    return print_expression(law)


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

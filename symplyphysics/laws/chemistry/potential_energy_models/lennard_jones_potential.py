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
## The Lennard-Jones (LJ) potential is an intermolecular pair potential. It models soft repulsive and
## attractive (van der Waals) interactions, and it describes electronically neutral atoms or molecules.
## It is a simple yet realistic model to describe intermolecular interactions: two particles repel each
## other at a very close distance, attract each other at moderate distance, and do not interact at infinite
## distance.

# Law: U = 4*e*((sigma/r)**12 - (sigma/r)**6)
## U - Lennard-Jones potential
## e - depth of potential well ("dispersion energy")
## sigma - distance at which the potential energy is zero ("particle size")
## r - distance between interacting particles

potential = Symbol("potential", units.energy)
dispersion_energy = Symbol("dispersion_energy", units.energy)
particle_size = Symbol("particle_size", units.length)
distance = Symbol("distance", units.length)

law = Eq(
    potential, 
    4 * dispersion_energy * ((particle_size / distance)**12 - (particle_size / distance)**6)
)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    dispersion_energy_=dispersion_energy,
    particle_size_=particle_size,
    distance_=distance,
)
@validate_output(potential)
def calculate_potential(
    dispersion_energy_: Quantity,
    particle_size_: Quantity,
    distance_: Quantity,
):
    result = law.rhs.subs({
        dispersion_energy: dispersion_energy_,
        particle_size: particle_size_,
        distance: distance_,
    })
    return Quantity(result)

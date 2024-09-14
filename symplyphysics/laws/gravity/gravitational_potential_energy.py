from sympy import Eq
from sympy.physics.units import gravitational_constant
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    clone_as_symbol,
    symbols,
)

# Description
## The gravitational potential energy U(r) of a system of two particles, with masses m1 and m2
## and separated by a distance r, is the negative of the work that would be done by the gravitational
## force of either particle acting on the other if the separation between the particles were changed
## from infinite to r.

# Law: U = -G * m1 * m2 / r
## U - gravitational potential energy
## G - gravitational constant
## m1, m2 - masses of interacting particles
## r - distance between mass centers of particles

gravitational_potential_energy = Symbol("gravitational_potential_energy", units.energy)
first_mass = clone_as_symbol(symbols.mass, display_symbol="m_1")
second_mass = clone_as_symbol(symbols.mass, display_symbol="m_2")
distance_between_mass_centers = Symbol("distance_between_mass_centers", units.length)

law = Eq(gravitational_potential_energy,
    -1 * gravitational_constant * first_mass * second_mass / distance_between_mass_centers)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    first_mass_=first_mass,
    second_mass_=second_mass,
    distance_between_mass_centers_=distance_between_mass_centers,
)
@validate_output(gravitational_potential_energy)
def calculate_gravitational_potential_energy(
    first_mass_: Quantity,
    second_mass_: Quantity,
    distance_between_mass_centers_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        first_mass: first_mass_,
        second_mass: second_mass_,
        distance_between_mass_centers: distance_between_mass_centers_,
    })
    return Quantity(result)

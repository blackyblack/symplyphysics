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
## Molar mass of a substance is the mass of 1 mole of particles that it is comprised of.

molar_mass = Symbol("molas_mass", units.mass / units.amount_of_substance)
particle_mass = Symbol("particle_mass", units.mass)

law = Eq(molar_mass, particle_mass * units.avogadro)


def print_law() -> str:
    return print_expression(law)


@validate_input(particle_mass_=particle_mass)
@validate_output(molar_mass)
def calculate_molar_mass(particle_mass_: Quantity) -> Quantity:
    result = law.rhs.subs(particle_mass, particle_mass_)
    return Quantity(result)

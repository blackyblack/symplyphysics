from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    symbols,
    clone_symbol,
)

# Description
## Rotational inertia of a rotating particle is defined as the product of its mass
## and the squared radius of rotation.

# Law: I = m * r^2
## I - rotational inertia
## m - particle mass
## r - radius of rotation

rotational_inertia = Symbol("rotational_inertia", units.mass * units.length**2)
radius = Symbol("radius", units.length)
particle_mass = clone_symbol(symbols.basic.mass, "particle_mass")

law = Eq(rotational_inertia, particle_mass * radius**2)


def print_law() -> str:
    return print_expression(law)


@validate_input(mass_=particle_mass, radius_=radius)
@validate_output(rotational_inertia)
def calculate_rotational_inertia(mass_: Quantity, radius_: Quantity) -> Quantity:
    result = law.rhs.subs({
        particle_mass: mass_,
        radius: radius_,
    })
    return Quantity(result)

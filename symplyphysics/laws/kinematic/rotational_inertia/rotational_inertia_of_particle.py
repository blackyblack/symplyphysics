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
## Rotational inertia of a rotating particle is defined as the product of its mass
## and the squared radius of rotation.

# Law: I = m * r^2
## I - rotational inertia
## m - particle mass
## r - radius of rotation

rotational_inertia = Symbol("rotational_inertia", units.mass * units.length**2)
mass = Symbol("mass", units.mass)
radius = Symbol("radius", units.length)

law = Eq(rotational_inertia, mass * radius**2)


def print_law() -> str:
    return print_expression(law)


@validate_input(mass_=mass, radius_=radius)
@validate_output(rotational_inertia)
def calculate_rotational_inertia(mass_: Quantity, radius_: Quantity) -> Quantity:
    result = law.rhs.subs({
        mass: mass_,
        radius: radius_,
    })
    return Quantity(result)

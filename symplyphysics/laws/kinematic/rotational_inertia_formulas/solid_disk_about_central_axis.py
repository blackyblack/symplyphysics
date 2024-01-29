from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)
from symplyphysics.laws.kinematic import rotational_inertia_cylindrical_integral as integral_law

# Description
## A solid disk (cylinder) rotates about its central axis (axis of cylindrical symmetry).

# Law: I = 1/2 * M * R**2
## I - rotational inertia
## M - mass of disk
## R - radius of disk

rotational_inertia = Symbol("rotational_inertia", units.mass * units.length**2)
mass = Symbol("mass", units.mass)
radius = Symbol("radius", units.length)

law = Eq(rotational_inertia, mass * radius**2 / 2)


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

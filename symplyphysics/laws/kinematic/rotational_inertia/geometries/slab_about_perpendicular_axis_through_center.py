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
## A solid slab is rotating about an axis that is perpendicular to the slab and goes through its center.

# Law: I = 1/12 * M * (a**2 + b**2)
## I - rotational inertia of slab
## M - mass of slab
## a, b - sizes of slab in the directions perpendicular to rotation axis

rotational_inertia = Symbol("rotational_inertia", units.mass * units.length**2)
mass = Symbol("mass", units.mass)
length = Symbol("length", units.length)
width = Symbol("width", units.length)

law = Eq(rotational_inertia, mass * (length**2 + width**2) / 12)


def print_law() -> str:
    return print_expression(law)


@validate_input(mass_=mass, length_=length, width_=width)
@validate_output(rotational_inertia)
def calculate_rotational_inertia(mass_: Quantity, length_: Quantity, width_: Quantity) -> Quantity:
    result = law.rhs.subs({
        mass: mass_,
        length: length_,
        width: width_,
    })
    return Quantity(result)

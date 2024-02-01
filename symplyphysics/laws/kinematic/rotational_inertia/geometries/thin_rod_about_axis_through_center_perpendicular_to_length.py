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
## A thin rod is rotating about an axis that passes through its center perpendicular to its length.

# Law: I = 1/12 * M * L**2
## I - rotational inertia of thin rod
## M - mass of rod
## L - length of rod

rotational_inertia = Symbol("rotational_inertia", units.mass * units.length**2)
mass = Symbol("mass", units.mass)
length = Symbol("length", units.length)

law = Eq(rotational_inertia, mass * length**2 / 12)


def print_law() -> str:
    return print_expression(law)


@validate_input(mass_=mass, length_=length)
@validate_output(rotational_inertia)
def calculate_rotational_inertia(mass_: Quantity, length_: Quantity) -> Quantity:
    result = law.rhs.subs({
        mass: mass_,
        length: length_,
    })
    return Quantity(result)

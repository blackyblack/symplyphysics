from sympy import Eq
from symplyphysics import (
    units,
    Symbol,
    Quantity,
    validate_input,
    validate_output,
    print_expression,
)

# Description
## The parallel-axis theorem relates the rotational inertia of a body about any axis to
## that of the same body about a parallel axis that extends through the body's center of mass
## of mass).

# Law: I = I_com + M*(h**2)
## I - rotational inertia about some axis
## I_com - rotational inertia about such an axis that is parallel to the given one
##         and that passes through the center of mass
## M - body mass
## h - perpendicular distance between the two axes

# Conditions:
## - The two axis must be parallel to each other.
## - The axis used in the calculation of I_com must pass through the body's center of mass.

rotational_inertia = Symbol("rotational_inertia", units.mass * units.length**2)
rotational_inertia_through_com = Symbol("rotational_inertia_through_com",
    units.mass * units.length**2)
mass = Symbol("mass", units.mass)
distance_between_axes = Symbol("distance_between_axes", units.length)

law = Eq(
    rotational_inertia,
    rotational_inertia_through_com + mass * distance_between_axes**2,
)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    rotational_inertia_through_com_=rotational_inertia_through_com,
    mass_=mass,
    distance_between_axes_=distance_between_axes,
)
@validate_output(rotational_inertia)
def calculate_rotational_inertia(
    rotational_inertia_through_com_: Quantity,
    mass_: Quantity,
    distance_between_axes_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        rotational_inertia_through_com: rotational_inertia_through_com_,
        mass: mass_,
        distance_between_axes: distance_between_axes_,
    })
    return Quantity(result)

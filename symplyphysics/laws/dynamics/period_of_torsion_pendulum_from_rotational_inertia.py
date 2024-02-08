from sympy import Eq, pi, sqrt
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## A torsion pendulum is an angular version of a linear harmonic oscillator: a disk
## oscillates in a horizontal plane; the reference line oscillates with some angular amplitude.
## The element of elasticity is associated with the twisting of the suspension wire.

# Law: T = 2 * pi * sqrt(I / kappa)
## T - period of the torsion pendulum
## I - rotational inertia of the disk
## kappa - torsion constant, depends on the properties of the suspension wire

oscillation_period = Symbol("oscillation_period", units.time)
rotational_inertia = Symbol("rotational_inertia", units.mass * units.length**2)
torsion_constant = Symbol("torsion_constant", units.force * units.length)

law = Eq(oscillation_period, 2 * pi * sqrt(rotational_inertia / torsion_constant))


# TODO: derive from relation between restoring torque and twist angle


def print_law() -> str:
    return print_expression(law)


@validate_input(rotational_inertia_=rotational_inertia, torsion_constant_=torsion_constant)
@validate_output(oscillation_period)
def calculate_period(rotational_inertia_: Quantity, torsion_constant_: Quantity) -> Quantity:
    result = law.rhs.subs({
        rotational_inertia: rotational_inertia_,
        torsion_constant: torsion_constant_,
    })
    return Quantity(result)

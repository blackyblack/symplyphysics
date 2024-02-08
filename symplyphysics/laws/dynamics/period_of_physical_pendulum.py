from sympy import Eq, pi, sqrt
from sympy.physics.units import acceleration_due_to_gravity
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## A physical pendulum is a pendulum with an arbitrary distribution of mass that
## oscillates about a given pivot point.

# Law: T = 2 * pi * sqrt(I / (m * g * h))
## T - period of physical pendulum
## I - rotational inertia of pendulum
## m - mass of pendulum
## g - acceleration due to gravity
## h - distance between pivot and pendulum's center of mass

oscillation_period = Symbol("oscillation_period", units.time)
rotational_inertia = Symbol("rotational_inertia", units.mass * units.length**2)
pendulum_mass = Symbol("pendulum_mass", units.mass)
distance_to_pivot = Symbol("distance_to_pivot", units.length)

law = Eq(
    oscillation_period, 
    2 * pi * sqrt(rotational_inertia / (pendulum_mass * acceleration_due_to_gravity * distance_to_pivot))
)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    rotational_inertia_=rotational_inertia,
    pendulum_mass_=pendulum_mass,
    distance_to_pivot_=distance_to_pivot,
)
@validate_output(oscillation_period)
def calculate_period(
    rotational_inertia_: Quantity,
    pendulum_mass_: Quantity,
    distance_to_pivot_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        rotational_inertia: rotational_inertia_,
        pendulum_mass: pendulum_mass_,
        distance_to_pivot: distance_to_pivot_,
    })
    return Quantity(result)

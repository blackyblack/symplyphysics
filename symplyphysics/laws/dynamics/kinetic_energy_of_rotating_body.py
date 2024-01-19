from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    angle_type,
    validate_input,
    validate_output,
)

# Description
## The kinetic energy of a rotating rigid body is expressed as half the product of
## its rotational inertia and the square of its angular velocity. Compare it with the
## kinetic energy of a moving body, where rotational inertia is replaced with gravitational
## mass of the body and angular velocity with linear velocity.

# Law: K = I * w^2 / 2
## K - kinetic energy
## I - rotational inertia
## w - angular velocity

kinetic_energy = Symbol("kinetic_energy", units.energy)
rotational_inertia = Symbol("rotational_inertia", units.mass * units.length**2)
angular_velocity = Symbol("angular_velocity", angle_type / units.time)

law = Eq(kinetic_energy, rotational_inertia * angular_velocity**2 / 2)


def print_law() -> str:
    return print_expression(law)


@validate_input(rotational_inertia_=rotational_inertia, angular_velocity_=angular_velocity)
@validate_output(kinetic_energy)
def calculate_kinetic_energy(rotational_inertia_: Quantity, angular_velocity_: Quantity) -> Quantity:
    result = solve(law, kinetic_energy)[0]
    result_energy = result.subs({
        rotational_inertia: rotational_inertia_,
        angular_velocity: angular_velocity_,
    })
    return Quantity(result_energy)

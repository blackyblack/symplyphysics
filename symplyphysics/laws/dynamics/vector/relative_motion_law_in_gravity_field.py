from sympy import Expr
from symplyphysics import (
    Vector,
    QuantityVector,
    units,
    Quantity,
    validate_input,
    validate_output,
    symbols,
    scale_vector,
    add_cartesian_vectors,
    subtract_cartesian_vectors,
    vector_magnitude,
)

# Description
## TODO

# Law: a = g - a_cor + F / m
## a - vector of acceleration relative to moving frame
## g - vector of acceleration due to gravity
## a_cor - vector of Coriolis acceleration
## F - vector sum of all forces acting on the body, excluding the Earth's gravitational pull
## m - mass of body

mass = symbols.basic.mass


def relative_acceleration_law(
    acceleration_due_to_gravity_: Vector,
    coriolis_acceleration_: Vector,
    total_force_: Vector,
) -> Vector:
    return add_cartesian_vectors(
        acceleration_due_to_gravity_,
        subtract_cartesian_vectors(
            scale_vector(1 / mass, total_force_),
            coriolis_acceleration_,
        )
    )


# F = m * (a - g + a_cor)
def force_law(
    relative_acceleration_: Vector,
    acceleration_due_to_gravity_: Vector,
    coriolis_acceleration_: Vector,
) -> Vector:
    return scale_vector(
        mass,
        add_cartesian_vectors(
            relative_acceleration_,
            subtract_cartesian_vectors(
                coriolis_acceleration_,
                acceleration_due_to_gravity_,
            )
        ),
    )


# m = norm(F) / norm(a - g - 2 * cross(v, w))
def mass_law(
    relative_acceleration_: Vector,
    acceleration_due_to_gravity_: Vector,
    coriolis_acceleration_: Vector,
    total_force_: Vector,
) -> Expr:
    acceleration_ = force_law(
        relative_acceleration_,
        acceleration_due_to_gravity_,
        coriolis_acceleration_,
    )
    return vector_magnitude(total_force_) / vector_magnitude(acceleration_)


@validate_input(
    acceleration_due_to_gravity_=units.acceleration,
    coriolis_acceleration_=units.acceleration,
    total_force_=units.force,
    mass_=mass,
)
@validate_output(units.acceleration)
def calculate_relative_acceleration(
    acceleration_due_to_gravity_: QuantityVector,
    coriolis_acceleration_: QuantityVector,
    total_force_: QuantityVector,
    mass_: Quantity,
) -> QuantityVector:
    vector = relative_acceleration_law(
        acceleration_due_to_gravity_.to_base_vector(),
        coriolis_acceleration_.to_base_vector(),
        total_force_.to_base_vector(),
    )
    return QuantityVector.from_base_vector(vector, subs={mass: mass_})

from sympy import Expr
from symplyphysics import (
    symbols,
    units,
    angle_type,
    Quantity,
    Vector,
    QuantityVector,
    validate_input,
    validate_output,
    scale_vector,
    add_cartesian_vectors,
    subtract_cartesian_vectors,
    cross_cartesian_vectors,
    vector_magnitude,
)

# Description
## The acceleration due to gravity is the acceleration a body in the gravity field of another body (e.g. Earth)
## with rotational effects such as the centripetal acceleration taken into account. It is the same for all bodies
## at a fixed point, but can have different magnitude and direction at different points in space.

# Law: g = F_gravity / m - cross(w, cross(w, r))
## g - vector of acceleration due to gravity
## F_gravity - vector of the force of gravity pull
## w - pseudovector of angular velocity of rotation of moving frame S'
## r - position vector in moving frame S'
## m - mass of body

mass = symbols.basic.mass


def acceleraton_due_to_gravity_law(
    gravity_force_: Vector,
    angular_velocity_: Vector,
    position_: Vector,
) -> Vector:
    return subtract_cartesian_vectors(
        scale_vector(1 / mass, gravity_force_),
        cross_cartesian_vectors(
            angular_velocity_,
            cross_cartesian_vectors(angular_velocity_, position_),
        ),
    )


def gravity_force_law(
    acceleration_due_to_gravity_: Vector,
    angular_velocity_: Vector,
    position_: Vector,
) -> Vector:
    return scale_vector(
        mass,
        add_cartesian_vectors(
            acceleration_due_to_gravity_,
            cross_cartesian_vectors(
                angular_velocity_,
                cross_cartesian_vectors(angular_velocity_, position_),
            ),
        ),
    )


def mass_law(
    acceleration_due_to_gravity_: Vector,
    gravity_force_: Vector,
    angular_velocity_: Vector,
    position_: Vector,
) -> Expr:
    acceleration_ = add_cartesian_vectors(
        acceleration_due_to_gravity_,
        cross_cartesian_vectors(
            angular_velocity_,
            cross_cartesian_vectors(angular_velocity_, position_),
        ),
    )
    return vector_magnitude(gravity_force_) / vector_magnitude(acceleration_)


@validate_input(
    gravity_force_=units.force,
    angular_velocity_=angle_type / units.time,
    position_=units.length,
    mass_=mass,
)
@validate_output(units.acceleration)
def calculate_acceleraton_due_to_gravity(
    gravity_force_: QuantityVector,
    angular_velocity_: QuantityVector,
    position_: QuantityVector,
    mass_: Quantity,
) -> QuantityVector:
    vector = acceleraton_due_to_gravity_law(
        gravity_force_.to_base_vector(),
        angular_velocity_.to_base_vector(),
        position_.to_base_vector(),
    )
    return QuantityVector.from_base_vector(vector, subs={mass: mass_})

from sympy import Expr
from symplyphysics import (
    units,
    symbols,
    Vector,
    Quantity,
    QuantityVector,
    validate_input,
    validate_output,
    scale_vector,
    add_cartesian_vectors,
    subtract_cartesian_vectors,
    vector_magnitude,
)

# Description
## TODO

# Law: a = g - a_cor + F / m
## a - vector of acceleration of body B in S'
## g - vector of acceleration due to gravity of body B
## a_cor - vector of Coriolis acceleration of body B in S'
## F - vector sum of all non-gravitational forces acting on body B
## m - mass of body B

mass = symbols.basic.mass


def acceleration_law(
    acceleration_due_to_gravity_: Vector,
    coriolis_acceleration_: Vector,
    force_: Vector,
) -> Vector:
    return subtract_cartesian_vectors(
        add_cartesian_vectors(
            acceleration_due_to_gravity_,
            scale_vector(1 / mass, force_),
        ),
        coriolis_acceleration_,
    )


# g = a + a_cor - F / m
def acceleration_due_to_gravity_law(
    acceleration_: Vector,
    coriolis_acceleration_: Vector,
    force_: Vector,
) -> Vector:
    return subtract_cartesian_vectors(
        add_cartesian_vectors(
            acceleration_,
            coriolis_acceleration_,
        ),
        scale_vector(1 / mass, force_),
    )


# a_cor = g - a + F / m
def coriolis_acceleration_law(
    acceleration_: Vector,
    acceleration_due_to_gravity_: Vector,
    force_: Vector,
) -> Vector:
    return subtract_cartesian_vectors(
        add_cartesian_vectors(
            acceleration_due_to_gravity_,
            scale_vector(1 / mass, force_),
        ),
        acceleration_,
    )


# F = m * (a + a_cor - g)
def force_law(
    acceleration_: Vector,
    acceleration_due_to_gravity_: Vector,
    coriolis_acceleration_: Vector,
) -> Vector:
    return scale_vector(
        mass,
        subtract_cartesian_vectors(
            add_cartesian_vectors(
                acceleration_,
                coriolis_acceleration_,
            ),
            acceleration_due_to_gravity_,
        ),
    )


# m = norm(F) / norm(a + a_cor - g) where `norm(x)` is Euclidean norm of vector `x`
def mass_law(
    acceleration_: Vector,
    acceleration_due_to_gravity_: Vector,
    coriolis_acceleration_: Vector,
    force_: Vector,
) -> Expr:
    total_acceleration_ = subtract_cartesian_vectors(
        add_cartesian_vectors(
            acceleration_,
            coriolis_acceleration_,
        ),
        acceleration_due_to_gravity_,
    )
    return vector_magnitude(force_) / vector_magnitude(total_acceleration_)


@validate_input(
    acceleration_due_to_gravity_=units.acceleration,
    coriolis_acceleration_=units.acceleration,
    force_=units.force,
    mass_=mass,
)
@validate_output(units.acceleration)
def calculate_acceleration(
    acceleration_due_to_gravity_: QuantityVector,
    coriolis_acceleration_: QuantityVector,
    force_: QuantityVector,
    mass_: Quantity,
) -> QuantityVector:
    vector = acceleration_law(
        acceleration_due_to_gravity_.to_base_vector(),
        coriolis_acceleration_.to_base_vector(),
        force_.to_base_vector(),
    )
    return QuantityVector.from_base_vector(
        vector,
        subs={mass: mass_},
    )

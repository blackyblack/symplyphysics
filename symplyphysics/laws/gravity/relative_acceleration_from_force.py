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
## Suppose reference frame S' is fixed to a moving body A (e.g. Earth). For some body B we can write a vector
## equation of motion relative to S' in the gravitational field of body A with body A's rotation taken into
## consideration. From this, we can gather the meaning of the acceleration due to gravity, also known as the
## free fall acceleration: it is the acceleration of body B relative to S' in the absence of external forces
## (`F = 0`) in the stationary case (the velocity of body B relative to S' is zero, i.e. `v = 0` and `a_cor = 0`).

# Law: a = g - a_cor + F / m
## a - vector of relative acceleration of body B in S'
## g - vector of acceleration of body B due due to the gravitational field of body A
## a_cor - vector of Coriolis acceleration of body B in S'
## F - vector sum of all non-gravitational forces acting on body B
## m - mass of body B

mass = symbols.mass


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

from sympy import Expr
from symplyphysics import (
    units,
    symbols,
    Quantity,
    Vector,
    QuantityVector,
    validate_input,
    validate_output,
    scale_vector,
    subtract_cartesian_vectors,
    add_cartesian_vectors,
    vector_magnitude,
)

# Description
## Suppose reference S' is fixed to a moving object (e.g. Earth). For some body B we can write an equation of motion
## in coordinates of S' akin to the Newton's second law of motion for inertial frames, although we obtain two
## additional components to the equation: one corresponding to the Coriolis force, and anther to the fictitious force
## of translation between inertial frame (S) and non-inertial frame S'.

# Law: a_rel = F / m - a_cor - a_tr
## a_rel - vector of acceleration of body B relative to S'
## F - vector sum of forces acting on body B
## a_cor - vector of [Coriolis acceleration](../../kinematic/vector/coriolis_acceleration.py) of body B
## a_tr - vector of [translation acceleration](../../kinematic/vector/acceleration_of_transfer_between_relative_frames.py) of body B

mass = symbols.basic.mass


def relative_acceleration_law(
    force_: Vector,
    coriolis_acceleration_: Vector,
    translation_acceleration_: Vector,
) -> Vector:
    return subtract_cartesian_vectors(
        subtract_cartesian_vectors(
            scale_vector(1 / mass, force_),
            coriolis_acceleration_,
        ),
        translation_acceleration_,
    )


# F = m * (a_rel + a_cor + a_tr)
def force_law(
    relative_acceleration_: Vector,
    coriolis_acceleration_: Vector,
    translation_acceleration_: Vector,
) -> Vector:
    return scale_vector(
        mass,
        add_cartesian_vectors(
            add_cartesian_vectors(
                relative_acceleration_,
                coriolis_acceleration_,
            ),
            translation_acceleration_,
        ),
    )


# m = norm(F) / norm(a_rel + a_cor + a_tr), where `norm(x)` is Euclidean norm of vector `x`
def mass_law(
    relative_acceleration_: Vector,
    force_: Vector,
    coriolis_acceleration_: Vector,
    translation_acceleration_: Vector,
) -> Expr:
    acceleration_ = add_cartesian_vectors(
        add_cartesian_vectors(
            relative_acceleration_,
            coriolis_acceleration_,
        ),
        translation_acceleration_,
    )
    return vector_magnitude(force_) / vector_magnitude(acceleration_)


# a_cor = F / m - a_rel - a_tr
def coriolis_acceleration_law(
    relative_acceleration_: Vector,
    force_: Vector,
    translation_acceleration_: Vector,
) -> Vector:
    return subtract_cartesian_vectors(
        subtract_cartesian_vectors(
            scale_vector(1 / mass, force_),
            relative_acceleration_,
        ),
        translation_acceleration_,
    )


# a_tr = F / m - a_rel - a_cor
def traslation_acceleration_law(
    relative_acceleration_: Vector,
    force_: Vector,
    coriolis_acceleration_: Vector,
) -> Vector:
    return subtract_cartesian_vectors(
        subtract_cartesian_vectors(
            scale_vector(1 / mass, force_),
            relative_acceleration_,
        ),
        coriolis_acceleration_,
    )


@validate_input(
    mass_=mass,
    force_=units.force,
    coriolis_acceleration_=units.acceleration,
    translation_acceleration_=units.acceleration,
)
@validate_output(units.acceleration)
def calculate_relative_acceleration(
    mass_: Quantity,
    force_: QuantityVector,
    coriolis_acceleration_: QuantityVector,
    translation_acceleration_: QuantityVector,
) -> QuantityVector:
    result_vector = relative_acceleration_law(
        force_.to_base_vector(),
        coriolis_acceleration_.to_base_vector(),
        translation_acceleration_.to_base_vector(),
    )
    result = QuantityVector.from_base_vector(
        result_vector,
        subs={mass: mass_},
    )
    return result

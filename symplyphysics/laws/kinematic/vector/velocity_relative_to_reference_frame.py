from sympy import Expr, symbols
from symplyphysics import (
    units,
    validate_input,
    validate_output,
    Vector,
    Quantity,
    QuantityVector,
    add_cartesian_vectors,
    subtract_cartesian_vectors,
    scale_vector,
)
from symplyphysics.core.vectors.arithmetics import (
    diff_cartesian_vector,
    integrate_cartesian_vector,
)

# Description
## For any reference frame, whether it is inertial or not, the motion relative to it can be described using
## the position vector relative to that frame's origin.

# Law: v_rel = (dx/dt, dy/dt, dz/dt)
## v_rel - vector of velocity relative to reference frame S
## r = (x, y, z) - vector of position relative to frame S
## d/dt - time derivative


def relative_velocity_law(
    position_: Vector,
    time_: Expr,
) -> Vector:
    return diff_cartesian_vector(position_, time_)


# r = r0 + Integral(v(t), t)
def relative_position_law(
    initial_position_: Vector,
    velocity_: Vector,
    time_: Expr,
) -> Vector:
    return add_cartesian_vectors(
        initial_position_,
        integrate_cartesian_vector(velocity_, time_),
    )


@validate_input(
    position_before_=units.length,
    position_after_=units.length,
    time_change_=units.time,
)
@validate_output(units.velocity)
def calculate_relative_velocity(
    position_before_: QuantityVector,
    position_after_: QuantityVector,
    time_change_: Quantity,
) -> QuantityVector:
    time_ = symbols("time")
    position_ = scale_vector(
        time_ / time_change_,
        subtract_cartesian_vectors(
            position_after_.to_base_vector(),
            position_before_.to_base_vector(),
        ),
    )
    velocity_ = relative_velocity_law(position_, time_)
    return QuantityVector.from_base_vector(velocity_)

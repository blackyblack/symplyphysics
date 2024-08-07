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

# Law: a_rel = (d**2(x)/dt**2, d**2(y)/dt**2, d**2(z)/dt**2)
## a_rel - vector of acceleration relative to reference frame S
## r = (x, y, z) - vector of position relative to frame S
## d**2/dt**2 - second-order derivative w.r.t. time

def relative_acceleration_law(
    position_: Vector,
    time_: Expr,
) -> Vector:
    return diff_cartesian_vector(position_, time_, time_)


# v = v0 + Integral(a(t), t)
def relative_velocity_law(
    initial_velocity_: Vector,
    acceleration_: Vector,
    time_: Expr,
) -> Vector:
    return add_cartesian_vectors(
        initial_velocity_,
        integrate_cartesian_vector(acceleration_, time_),
    )


# r = r0 + v0 * t + Integral(a(t), t, t)
def relative_position_law(
    initial_position_: Vector,
    initial_velocity_: Vector,
    acceleration_: Vector,
    time_: Expr,
) -> Vector:
    velocity_ = relative_velocity_law(initial_velocity_, acceleration_, time_)
    position_ = add_cartesian_vectors(
        initial_position_,
        integrate_cartesian_vector(velocity_, time_)
    )
    return position_


# (dr**2/dt**2)(t) ~= (r(t - dt) - 2 * r(t) + r(t + dt)) / dt**2
@validate_input(
    position_before_=units.length,
    position_current_=units.length,
    position_after_=units.length,
    time_change_=units.time,
)
@validate_output(units.acceleration)
def calculate_relative_acceleration(
    position_before_: Vector,
    position_current_: Vector,
    position_after_: Vector,
    time_change_: Quantity,
) -> QuantityVector:
    time_ = symbols("time")
    position_ = scale_vector(
        (time_ / time_change_)**2 / 2,
        add_cartesian_vectors(
            subtract_cartesian_vectors(
                position_before_,
                scale_vector(2, position_current_),
            ),
            position_after_,
        )
    )
    acceleration_ = relative_acceleration_law(position_, time_)
    return QuantityVector.from_base_vector(acceleration_)

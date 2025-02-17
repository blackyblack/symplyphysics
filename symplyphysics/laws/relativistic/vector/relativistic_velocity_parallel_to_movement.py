from typing import Optional
from sympy.physics.units import speed_of_light
from symplyphysics import (
    units,
    Quantity,
    validate_input,
    validate_output,
    Vector,
    dot_vectors,
    add_cartesian_vectors,
    scale_vector,
    QuantityVector,
    cross_cartesian_vectors,
    vector_magnitude,
    assert_equal,
)

# Description
## Consider two inertial reference frames: one fixed (lab frame) and one tied to the moving object (proper frame).
## The proper frame is moving with some velocity `v` relative to the lab frame. Then, according to the theory of special
## relativity, the velocity of the object relative to lab frame is not equal to the sum of its velocity in the proper
## frame and the velocity of the proper frame relative to the lab frame.

# Law: u_parallel = (u_parallel' + v) / (1 + dot(u_parallel', v) / c**2)
## u_parallel - velocity vector relative to lab frame parallel to `v`
## u_parallel' - velocity vector relative to proper frame parallel to `v`
## v - velocity vector of proper frame relative to lab frame
## c - speed of light

# Notes
## - One can get the same expression for `u_parallel'` in terms of `u_parallel` by replacing `v` with `-v`. This is
##   essentially the inverse Lorentz transformation from lab frame to proper frame that uses the fact that the lab frame
##   can be viewed as moving with velocity vector `-v` relative to the proper frame.

# Conditions
## - Works in special relativity

# Links: Wikipedia <https://en.wikipedia.org/wiki/Velocity-addition_formula#General_configuration>


def parallel_velocity_component_in_lab_frame_law(
    parallel_velocity_component_in_proper_frame_: Vector,
    proper_frame_velocity_: Vector,
) -> Vector:
    dot = dot_vectors(
        parallel_velocity_component_in_proper_frame_,
        proper_frame_velocity_,
    )

    factor = 1 + dot / speed_of_light**2

    added = add_cartesian_vectors(
        parallel_velocity_component_in_proper_frame_,
        proper_frame_velocity_,
    )

    return scale_vector(1 / factor, added)


def parallel_velocity_component_in_proper_frame_law(
    parallel_velocity_component_in_lab_frame_: Vector,
    lab_frame_velocity_: Vector,
) -> Vector:
    return parallel_velocity_component_in_lab_frame_law(
        parallel_velocity_component_in_lab_frame_,
        scale_vector(-1, lab_frame_velocity_),
    )


def proper_frame_velocity_in_lab_frame_law(
    parallel_velocity_component_in_proper_frame_: Vector,
    parallel_velocity_component_in_lab_frame_: Vector,
) -> Vector:
    lab_dot_lab = dot_vectors(
        parallel_velocity_component_in_lab_frame_,
        parallel_velocity_component_in_lab_frame_,
    )

    proper_dot_lab = dot_vectors(
        parallel_velocity_component_in_proper_frame_,
        parallel_velocity_component_in_lab_frame_,
    )

    factor = (lab_dot_lab - proper_dot_lab) / (1 - proper_dot_lab / speed_of_light**2)

    return scale_vector(factor / lab_dot_lab, parallel_velocity_component_in_lab_frame_)


@validate_input(
    parallel_velocity_component_in_proper_frame_=units.velocity,
    proper_frame_velocity_=units.velocity,
)
@validate_output(units.velocity)
def calculate_parallel_velocity_component_in_lab_frame(
    parallel_velocity_component_in_proper_frame_: QuantityVector,
    proper_frame_velocity_: QuantityVector,
    *,
    tolerance_: Optional[float] = None,
) -> QuantityVector:
    cross = cross_cartesian_vectors(
        parallel_velocity_component_in_proper_frame_.to_base_vector(),
        proper_frame_velocity_.to_base_vector(),
    )

    assert_equal(Quantity(vector_magnitude(cross)).scale_factor, 0, relative_tolerance=tolerance_)

    result = parallel_velocity_component_in_lab_frame_law(
        parallel_velocity_component_in_proper_frame_.to_base_vector(),
        proper_frame_velocity_.to_base_vector(),
    )
    return QuantityVector.from_base_vector(result)

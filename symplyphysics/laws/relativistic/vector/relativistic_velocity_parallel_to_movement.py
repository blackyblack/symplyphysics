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
from symplyphysics.core.approx import APPROX_RELATIVE_TOLERANCE

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


@validate_input(
    parallel_velocity_component_in_proper_frame_=units.velocity,
    proper_frame_velocity_=units.velocity,
)
@validate_output(units.velocity)
def calculate_parallel_velocity_component_in_lab_frame(
    parallel_velocity_component_in_proper_frame_: QuantityVector,
    proper_frame_velocity_: QuantityVector,
    *,
    tolerance_: float = APPROX_RELATIVE_TOLERANCE,
) -> QuantityVector:
    cross = cross_cartesian_vectors(
        parallel_velocity_component_in_proper_frame_.to_base_vector(),
        proper_frame_velocity_.to_base_vector(),
    )

    assert_equal(Quantity(vector_magnitude(cross)).scale_factor, 0, tolerance=tolerance_)

    result = parallel_velocity_component_in_lab_frame_law(
        parallel_velocity_component_in_proper_frame_.to_base_vector(),
        proper_frame_velocity_.to_base_vector(),
    )
    return QuantityVector.from_base_vector(result)

from sympy.physics.units import speed_of_light
from symplyphysics import (
    units,
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
## TODO

# Law: u_parallel = (u_parallel' + v) / (1 + dot(u_parallel', v) / c**2)
## u_parallel - velocity vector relative to lab frame in the perpendicular direction
## u_parallel' - velocity vector relative to proper frame in the perpendicular direction
## v - velocity vector of proper frame relative to lab frame
## c - speed of light


def parallel_addition_law(
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

    scaled = scale_vector(factor, added)
    return scaled


@validate_input(
    parallel_velocity_component_in_proper_frame_=units.velocity,
    proper_frame_velocity_=units.velocity,
)
@validate_output(units.velocity)
def calculate_parallel_velocity_component_in_lab_frame(
    parallel_velocity_component_in_proper_frame_: QuantityVector,
    proper_frame_velocity_: QuantityVector,
) -> QuantityVector:
    cross = cross_cartesian_vectors(
        parallel_velocity_component_in_proper_frame_.to_base_vector(),
        proper_frame_velocity_.to_base_vector(),
    )

    try:
        assert_equal(vector_magnitude(cross), 0)
    except AssertionError as e:
        raise ValueError("The two velocity vectors must be collinear") from e

    result = parallel_addition_law(
        parallel_velocity_component_in_proper_frame_,
        proper_frame_velocity_,
    )
    return QuantityVector.from_base_vector(result)

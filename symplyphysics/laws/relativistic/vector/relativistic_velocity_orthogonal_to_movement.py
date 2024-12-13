from sympy.physics.units import speed_of_light
from symplyphysics import (
    units,
    validate_input,
    validate_output,
    Vector,
    vector_magnitude,
    dot_vectors,
    scale_vector,
    QuantityVector,
)
from symplyphysics.core.vectors.arithmetics import reject_cartesian_vector
from symplyphysics.definitions import lorentz_factor as lorentz_factor_def

# Description
## Consider two inertial reference frames: one fixed (lab frame) and one tied to the moving object (proper frame).
## The proper frame is moving with some velocity `v` relative to the lab frame. Then, according to the theory of special
## relativity, the velocity of the object relative to lab frame is not equal to the sum of its velocity in the proper
## frame and the velocity of the proper frame relative to the lab frame.

# Law: u_orthogonal = u_orthogonal' / (gamma * (1 + dot(v, u') / c**2))
## u_orthogonal - component of velocity vector relative to lab frame orthogonal to `v`
## u_orthogonal' - component of velocity vector relative to proper frame orthogonal to `v`
## u' - velocity vector relative to proper frame
## v - velocity vector of proper frame relative to lab frame
## gamma = 1 / sqrt(1 - dot(v, v) / c**2) - Lorentz factor
## dot(a, b) - dot product between vectors `a` and `b`

# Notes
## - One can get the same expression for `u_orthogonal'` in terms of `u_orthogonal` by replacing `v` with `-v`. This is
##   essentially the inverse Lorentz transformation from lab frame to proper frame that uses the fact that the lab frame
##   can be viewed as moving with velocity vector `-v` relative to the proper frame.

# Conditions
## - Works in special relativity

# Links: Wikipedia <https://en.wikipedia.org/wiki/Velocity-addition_formula#General_configuration>


def orthogonal_velocity_component_in_lab_frame_law(
    velocity_in_proper_frame_: Vector,
    proper_frame_velocity_: Vector,
) -> Vector:
    orthogonal_velocity_component_in_proper_frame_ = reject_cartesian_vector(
        velocity_in_proper_frame_, proper_frame_velocity_)

    lorentz_factor_ = lorentz_factor_def.definition.rhs.subs({
        lorentz_factor_def.speed: vector_magnitude(proper_frame_velocity_),
    })

    scale_factor_ = lorentz_factor_ * (1 +
        dot_vectors(proper_frame_velocity_, velocity_in_proper_frame_) / speed_of_light**2)

    return scale_vector(1 / scale_factor_, orthogonal_velocity_component_in_proper_frame_)


def orthogonal_velocity_component_in_proper_frame_law(
    velocity_in_lab_frame_: Vector,
    lab_frame_velocity_: Vector,
) -> Vector:
    return orthogonal_velocity_component_in_lab_frame_law(
        velocity_in_lab_frame_,
        scale_vector(-1, lab_frame_velocity_),
    )


@validate_input(
    velocity_in_proper_frame_=units.velocity,
    proper_frame_velocity_=units.velocity,
)
@validate_output(units.velocity)
def calculate_orthogonal_velocity_component_in_lab_frame(
    velocity_in_proper_frame_: QuantityVector,
    proper_frame_velocity_: QuantityVector,
) -> QuantityVector:
    result = orthogonal_velocity_component_in_lab_frame_law(
        velocity_in_proper_frame_.to_base_vector(),
        proper_frame_velocity_.to_base_vector(),
    )
    return QuantityVector.from_base_vector(result)

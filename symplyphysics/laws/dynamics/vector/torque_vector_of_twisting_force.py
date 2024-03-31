from symplyphysics import (
    units,
    QuantityVector,
    Vector,
    cross_cartesian_vectors,
    validate_input,
    validate_output,
)

# Description
## Torque is a turning or twisting action on a body about a rotation axis due to a force.
## It is a pseudovector defined as a cross product of the force vector and the position vector of the point
## of force application.

# Law: tau = cross(r, F)
## tau - torque pseudovector
## r - position vector of the point of force applicaton
## F - force vector
## cross(a, b) - cross product between vectors a and b


def torque_definition(position_: Vector, force_: Vector) -> Vector:
    return cross_cartesian_vectors(position_, force_)


@validate_input(force_=units.force, position_=units.length)
@validate_output(units.force * units.length)
def calculate_torque(position_: QuantityVector, force_: QuantityVector) -> QuantityVector:
    torque_vector = torque_definition(position_.to_base_vector(), force_.to_base_vector())
    return QuantityVector.from_base_vector(torque_vector)

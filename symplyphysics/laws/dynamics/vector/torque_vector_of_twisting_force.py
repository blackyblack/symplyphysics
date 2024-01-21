from symplyphysics import (
    units,
    QuantityVector,
    Vector,
    cross_cartesian_vectors,
    validate_input,
    validate_output,
    list_of_quantities,
)

# Description
## Torque is a turning or twisting action on a body about a rotation axis due to a force.
## It is a pseudovector defined as a cross product of the force vector and the position vector of the point
## of force application.

# Law: tau = [F, r]
## tau - torque vector
## F - force vector
## r - position vector of the point of force applicaton


def torque_definition(force_: Vector, position_: Vector) -> Vector:
    return cross_cartesian_vectors(force_, position_)


@validate_input(force_=units.force, position_=units.length)
@validate_output(units.force * units.length)
def calculate_torque(force_: QuantityVector, position_: QuantityVector) -> QuantityVector:
    torque_vector = torque_definition(force_, position_)
    return QuantityVector(torque_vector.components, force_.coordinate_system)

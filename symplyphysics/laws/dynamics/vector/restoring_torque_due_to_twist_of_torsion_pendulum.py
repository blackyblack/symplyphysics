from symplyphysics import (
    Vector,
    QuantityVector,
    scale_vector,
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    subs_list,
    angle_type,
)

# Description
## A torsion pendulum is a simple harmonic oscillator consisting of a disk suspended by a wire.
## Rotating the disk through an angle in either direction introduces a restoring torque.

# Law: tau = -kappa * theta
## tau - torque pseudovector
## kappa - torsion constant of the pendulum
## theta - angular displacement pseudovector, also known as rotation vector

torsion_constant = Symbol("torsion_constant", units.force * units.length)


def torque_law(rotation_vector_: Vector) -> Vector:
    return scale_vector(-1 * torsion_constant, rotation_vector_)


def rotation_vector_law(torque_: Vector) -> Vector:
    return scale_vector(-1 / torsion_constant, torque_)


@validate_input(torsion_constant_=torsion_constant, rotation_vector_=angle_type)
@validate_output(units.force * units.length)
def calculate_torque(torsion_constant_: Quantity,
    rotation_vector_: QuantityVector) -> QuantityVector:
    result_vector = torque_law(rotation_vector_.to_base_vector())
    result_components = subs_list(result_vector.components,
        {torsion_constant: torsion_constant_})
    return QuantityVector(result_components, rotation_vector_.coordinate_system)


@validate_input(torsion_constant_=torsion_constant, torque_=units.force * units.length)
@validate_output(angle_type)
def calculate_rotation_vector(torsion_constant_: Quantity,
    torque_: QuantityVector) -> QuantityVector:
    result_vector = rotation_vector_law(torque_.to_base_vector())
    result_components = subs_list(result_vector.components,
        {torsion_constant: torsion_constant_})
    return QuantityVector(result_components, torque_.coordinate_system)

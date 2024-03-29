from symplyphysics import (
    units,
    Symbol,
    Quantity,
    Vector,
    QuantityVector,
    scale_vector,
    validate_input,
    validate_output,
)

# Description
## Damping force is an external (relative to an object) force that drains energy from the object,
## reducing the motion of the object. It is a model used, for example, describing the motion
## of an oscillator.

# Law: F = -b * v
## F - vector of damping force
## b - damping constant, b >= 0
## v - vector of object's velocity

damping_constant = Symbol("damping_constant", units.mass / units.time)


def damping_force_definition(velocity_: Vector) -> Vector:
    return scale_vector(-1 * damping_constant, velocity_)


def velocity_law(damping_force_: Vector) -> Vector:
    return scale_vector(-1 / damping_constant, damping_force_)


@validate_input(damping_constant_=damping_constant, velocity_=units.velocity)
@validate_output(units.force)
def calculate_damping_force(damping_constant_: Quantity,
    velocity_: QuantityVector) -> QuantityVector:
    result_vector = damping_force_definition(velocity_.to_base_vector())
    return QuantityVector.from_base_vector(result_vector,
        subs={damping_constant: damping_constant_})


@validate_input(damping_constant_=damping_constant, damping_force_=units.force)
@validate_output(units.velocity)
def calculate_velocity(damping_constant_: Quantity,
    damping_force_: QuantityVector) -> QuantityVector:
    result_vector = velocity_law(damping_force_.to_base_vector())
    return QuantityVector.from_base_vector(result_vector,
        subs={damping_constant: damping_constant_})

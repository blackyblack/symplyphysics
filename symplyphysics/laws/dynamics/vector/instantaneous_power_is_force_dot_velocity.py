from sympy import Expr
from symplyphysics import units, Quantity, Vector, validate_input, validate_output
from symplyphysics.core.vectors.arithmetics import dot_vectors
from symplyphysics.core.vectors.vectors import QuantityVector

# Description
## The power due to a force is the rate at which that force does work on an object.

# Law: P = dot(F, v)
## P - instantaneous power
## F - force vector
## v - instantaneous velocity vector
## dot(a, b) is the dot product between vectors a and b


def power_law(force_: Vector, velocity_: Vector) -> Expr:
    return dot_vectors(force_, velocity_)


@validate_input(force_=units.force, velocity_=units.velocity)
@validate_output(units.power)
def calculate_power(force_: QuantityVector, velocity_: QuantityVector) -> Quantity:
    result = power_law(force_.to_base_vector(), velocity_.to_base_vector())
    return Quantity(result)

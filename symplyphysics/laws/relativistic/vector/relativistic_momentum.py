from sympy import sqrt
from sympy.physics.units import speed_of_light
from symplyphysics import (
    units,
    Symbol,
    Quantity,
    validate_input,
    validate_output,
    Vector,
    dot_vectors,
    scale_vector,
    QuantityVector,
)

# Description
## Momentum (amount of motion) is a vector physical quantity that is a measure of the mechanical movement
## of a body. The relativistic momentum also takes into account speed limits equal to the speed of light.

# Law: p = m0 * v / sqrt(1 - |v|**2 / c**2)
## p - vector of relativistic momentum
## m0 - rest (invariant) mass of object
## c - speed of light
## v - vector of particle's velocity
## |v| - magnitude of velocity

rest_mass = Symbol("rest_mass", units.mass)


def momentum_law(velocity_: Vector) -> Vector:
    v_dot_v = dot_vectors(velocity_, velocity_)
    factor = rest_mass / sqrt(1 - v_dot_v / speed_of_light**2)
    return scale_vector(factor, velocity_)


def velocity_law(momentum_: Vector) -> Vector:
    p_dot_p = dot_vectors(momentum_, momentum_)
    factor = speed_of_light / sqrt((rest_mass * speed_of_light)**2 + p_dot_p)
    return scale_vector(factor, momentum_)


@validate_input(
    rest_mass_=rest_mass,
    velocity_=units.velocity,
)
@validate_output(units.momentum)
def calculate_momentum(
    rest_mass_: Quantity,
    velocity_: QuantityVector,
) -> QuantityVector:
    result = momentum_law(velocity_.to_base_vector())
    return QuantityVector.from_base_vector(
        result,
        subs={rest_mass: rest_mass_},
    )

from sympy import Derivative
from symplyphysics import (
    units,
    Symbol,
    Quantity,
    Vector,
    validate_input,
    validate_output,
    scale_vector,
    add_cartesian_vectors,
)
from symplyphysics.core.vectors.vectors import QuantityVector

# Description
## Newton's second law of motion can be generalized in terms of the vector of linear momentum.

# Law:
## F = dp/dt
## F - force vector
## p - vector of linear momentum

time = Symbol("time", units.time)


def force_law(momentum_: Vector) -> Vector:
    force_components = list(map(
        lambda component: Derivative(component, time),
        momentum_.components,
    ))
    return Vector(force_components, momentum_.coordinate_system)


@validate_input(
    momentum_before_=units.mass * units.velocity,
    momentum_after_=units.mass * units.velocity,
    time_=time,
)
@validate_output(units.force)
def calculate_force(
    momentum_before_: QuantityVector,
    momentum_after_: QuantityVector,
    time_: Quantity,
) -> QuantityVector:
    # delta_p = (t / delta_t) * (p1 - p0)
    momentum_function = scale_vector(
        time / time_,
        add_cartesian_vectors(
            momentum_after_, 
            scale_vector(-1, momentum_before_)
        )
    )
    result = [component.doit() for component in force_law(momentum_function).components]  # type: ignore
    return QuantityVector(result, momentum_before_.coordinate_system)

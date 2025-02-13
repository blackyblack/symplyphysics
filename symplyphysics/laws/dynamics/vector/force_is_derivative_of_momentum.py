from sympy import Derivative, sympify
from symplyphysics import (
    units,
    Quantity,
    Vector,
    validate_input,
    validate_output,
    scale_vector,
    add_cartesian_vectors,
    symbols,
)
from symplyphysics.core.vectors.vectors import QuantityVector

# Description
## Newton's second law of motion can be generalized in terms of the vector of linear momentum.

# Law:
## F = dp/dt
## F - force vector
## p - vector of linear momentum

time = symbols.time


def force_law(momentum_: Vector) -> Vector:
    force_components = [Derivative(component, time) for component in momentum_.components]
    return Vector(force_components, momentum_.coordinate_system)


@validate_input(
    momentum_before_=units.momentum,
    momentum_after_=units.momentum,
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
        add_cartesian_vectors(momentum_after_.to_base_vector(),
        scale_vector(-1, momentum_before_.to_base_vector())))
    result_force_vector = force_law(momentum_function)
    result_force_components = [
        sympify(component).doit() for component in result_force_vector.components
    ]
    return QuantityVector(result_force_components, momentum_before_.coordinate_system)

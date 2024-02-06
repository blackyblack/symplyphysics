from sympy import Derivative, sympify
from symplyphysics import (
    units,
    Symbol,
    Quantity,
    QuantityVector,
    Vector,
    validate_input,
    validate_output,
    scale_vector,
    add_cartesian_vectors,
)

# Description
## - Case of a single particle
##   The total vector sum of all the external torques acting on a particle is equal to the 
##   time rate change of the angular momentum of that particle.
## - Case of a system of particles
##   The net external torque acting on a system of particles is equal to the time rate change
##   of the system's total angular momentum.

## Law: tau = dL/dt
## tau - vector of net torque
## L - vector of (total) angular momentum
## t - time

time = Symbol("time", units.time)


def torque_law(angular_momentum_: Vector) -> Vector:
    torque_components = list(map(
        lambda component: Derivative(component, time),
        angular_momentum_.components,
    ))
    return Vector(torque_components, angular_momentum_.coordinate_system)


@validate_input(
    angular_momentum_before_=units.length * units.momentum,
    angular_momentum_after_=units.length * units.momentum,
    time_=time,
)
@validate_output(units.force * units.length)
def calculate_torque(
    angular_momentum_before_: QuantityVector,
    angular_momentum_after_: QuantityVector,
    time_: Quantity,
) -> QuantityVector:
    angular_momentum_function = scale_vector(
        time / time_,
        add_cartesian_vectors(
            angular_momentum_after_,
            scale_vector(-1, angular_momentum_before_)
        )
    )
    result_components = [
        sympify(component).doit()
        for component in torque_law(angular_momentum_function).components
    ]
    return QuantityVector(result_components, angular_momentum_before_.coordinate_system)

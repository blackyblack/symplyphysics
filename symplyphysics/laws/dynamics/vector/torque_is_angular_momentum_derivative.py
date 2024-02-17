from sympy import Derivative, sympify
from symplyphysics import (
    CoordinateSystem,
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
    torque_components = [Derivative(component, time) for component in angular_momentum_.components]
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
    if angular_momentum_before_.coordinate_system.coord_system_type != CoordinateSystem.System.CARTESIAN:
        raise ValueError("Initial angular momentum vector should be in cartesian coordinate system")
    if angular_momentum_after_.coordinate_system.coord_system_type != CoordinateSystem.System.CARTESIAN:
        raise ValueError("Final angular momentum vector should be in cartesian coordinate system")
    angular_momentum_function = scale_vector(
        time / time_,
        add_cartesian_vectors(angular_momentum_after_.to_base_vector(),
        scale_vector(-1, angular_momentum_before_.to_base_vector())))
    result_torque_vector = torque_law(angular_momentum_function)
    result_components = [sympify(component).doit() for component in result_torque_vector.components]
    return QuantityVector(result_components, angular_momentum_before_.coordinate_system)

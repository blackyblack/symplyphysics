"""
Torque is angular momentum derivative
=====================================

In the case of a single particle, the total vector sum of all the external torques acting on a
particle is equal to the time rate change of the angular momentum of that particle.

In the case of a system of particles, the net external torque acting on a system of particles is
equal to the time rate change of the system's total angular momentum.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Torque#Relationship_with_the_angular_momentum>`__.
"""

from sympy import Eq
from symplyphysics import Quantity, validate_input, validate_output, symbols

from symplyphysics.core.experimental.vectors import (clone_as_vector_function, VectorDerivative,
    vector_diff, convert_sympy_to_vector_derivatives)
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector
from symplyphysics.core.experimental.solvers import solve_for_vector, apply, vector_equals

from symplyphysics.definitions.vector import (
    velocity_is_position_vector_derivative as _velocity_def,
    angular_momentum_is_position_cross_linear_momentum as _angular_momentum_def,
    momentum_is_mass_times_velocity_vector as _linear_momentum_def,
    acceleration_is_velocity_derivative as _acceleration_def,
)
from symplyphysics.laws.dynamics.vector import (
    acceleration_from_force as _newtons_second_law,
    torque_vector_of_twisting_force as _torque_def,
)

time = symbols.time
"""
:symbols:`time`.
"""

angular_momentum = clone_as_vector_function(symbols.angular_momentum, (time,))
"""
Pseudovector of net :symbols:`angular_momentum` as a function of :attr:`~time`.
"""

torque = clone_as_vector_function(symbols.torque, (time,))
"""
Pseudovector of net :symbols:`torque` as a function of :attr:`~time`.
"""

law = Eq(torque(time), VectorDerivative(angular_momentum(time), time))
"""
:laws:symbol::

:laws:latex::
"""

# Derive law

_mass = _newtons_second_law.mass

_position = _velocity_def.position_vector

_velocity_eqn = _velocity_def.law.subs(_velocity_def.time, time).subs(
    _velocity_def.position_vector(time),
    _position(time),
)

_velocity_expr = solve_for_vector(_velocity_eqn, _velocity_def.velocity(time))

_linear_momentum_eqn = _linear_momentum_def.law.subs({
    _linear_momentum_def.mass: _mass,
    _linear_momentum_def.velocity: _velocity_expr,
})

_linear_momentum_expr = solve_for_vector(_linear_momentum_eqn, _linear_momentum_def.momentum)

_angular_momentum_eqn = _angular_momentum_def.law.subs({
    _angular_momentum_def.angular_momentum: angular_momentum(time),
    _angular_momentum_def.position_vector: _position(time),
    _angular_momentum_def.linear_momentum: _linear_momentum_expr,
})

_acceleration_eqn = _acceleration_def.law.subs(_acceleration_def.time, time).subs({
    _acceleration_def.acceleration(time): _newtons_second_law.acceleration,
    _acceleration_def.velocity(time): _velocity_expr,
})
_acceleration_eqn = apply(_acceleration_eqn, convert_sympy_to_vector_derivatives)

_newtons_second_eqn = _newtons_second_law.law.subs(
    _newtons_second_law.force,
    _torque_def.force,
)

_torque_eqn = _torque_def.law.subs({
    _torque_def.torque: torque(time),
    _torque_def.position_vector: _position(time),
})

_derived_law = apply(
    _angular_momentum_eqn,
    lambda x: vector_diff(x, time),
).subs(
    _acceleration_eqn.rhs,
    _acceleration_eqn.lhs,
).subs(
    _newtons_second_eqn.lhs,
    _newtons_second_eqn.rhs,
).subs(
    _torque_eqn.rhs,
    _torque_eqn.lhs,
)

# The sides are swapped in `_derived_law`
assert vector_equals(_derived_law.rhs, law.lhs)
assert vector_equals(_derived_law.lhs, law.rhs)


@validate_input(
    angular_momentum_before_=angular_momentum,
    angular_momentum_after_=angular_momentum,
    time_=time,
)
@validate_output(torque)
def calculate_torque(
    angular_momentum_before_: QuantityCoordinateVector,
    angular_momentum_after_: QuantityCoordinateVector,
    time_: Quantity,
) -> QuantityCoordinateVector:
    angular_momentum_ = (time / time_) * (angular_momentum_after_ - angular_momentum_before_)

    result = solve_for_vector(
        law.subs(angular_momentum(time), angular_momentum_),
        torque(time),
    ).doit()

    return QuantityCoordinateVector.from_expr(result)

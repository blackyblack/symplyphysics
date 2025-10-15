"""
Instantaneous power is dot product of force and velocity
========================================================

**Instantaneous power** is a measure of the rate of energy transfer or conversion at a given point
in time.

**Conditions:**

#. Force :math:`\\vec F` is constant.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Power_(physics)#Definition>`__.
"""

from sympy import Eq

from symplyphysics import Quantity, validate_input, validate_output, symbols
from symplyphysics.core.expr_comparisons import expr_equals

from symplyphysics.core.experimental.vectors import (clone_as_vector_symbol, VectorDot,
    clone_as_vector_function, convert_sympy_to_vector_derivatives)
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

from symplyphysics.definitions import power_is_energy_derivative as _power_def
from symplyphysics.definitions.vector import velocity_is_position_vector_derivative as _velocity_def
from symplyphysics.laws.dynamics.vector import mechanical_work_from_force_and_displacement as _work_law

power = symbols.power
"""
Instantaneous :symbols:`power` due to :attr:`~force`.
"""

force = clone_as_vector_symbol(symbols.force)
"""
Vector of :symbols:`force` exerted on the body.
"""

velocity = clone_as_vector_symbol(symbols.speed)
"""
Vector of the body's velocity. See :symbols:`speed`.
"""

law = Eq(power, VectorDot(force, velocity))
"""
:laws:symbol::

:laws:latex::
"""

# Derive law

_time = _power_def.time

_displacement = clone_as_vector_function(symbols.distance, (_time,))

_work_expr = _work_law.law.rhs.subs({
    _work_law.force: force,
    _work_law.displacement: _displacement(_time),
}).doit()

_power_expr = _power_def.definition.rhs.subs(
    _power_def.energy(_time),
    _work_expr,
).doit()

_velocity_eqn = _velocity_def.law.subs(_velocity_def.time, _time).subs({
    _velocity_def.position_vector(_time): _displacement(_time),
    _velocity_def.velocity(_time): velocity,
})

_power_expr = _power_expr.subs(
    convert_sympy_to_vector_derivatives(_velocity_eqn.rhs),
    _velocity_eqn.lhs,
)

assert expr_equals(_power_expr, law.rhs)


@validate_input(force_=force, velocity_=velocity)
@validate_output(power)
def calculate_power(
    force_: QuantityCoordinateVector,
    velocity_: QuantityCoordinateVector,
) -> Quantity:
    result = law.rhs.subs({
        force: force_,
        velocity: velocity_,
    })
    return Quantity(result)


# UNIQUE_LAW_ID: 239

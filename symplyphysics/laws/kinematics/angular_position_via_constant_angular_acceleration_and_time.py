r"""
Angular position via constant angular acceleration and time
===========================================================

If a body is rotating with a constant acceleration, its angular position is a quadratic function of time.

#. The axis is fixed.
#. Angular acceleration is constant, i.e. :math:`\frac{d \alpha}{d t} = 0.`
"""

from sympy import Eq, solve, dsolve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.symbols.quantities import scale_factor
from symplyphysics.definitions import (
    angular_acceleration_is_angular_speed_derivative as angular_acceleration_def,
    angular_speed_is_angular_distance_derivative as angular_velocity_def,
)

final_angular_position = symbols.angular_distance
"""
:symbols:`angular_distance` at :attr:`~time`.
"""

initial_angular_position = clone_as_symbol(symbols.angular_distance, subscript="0")
"""
:symbols:`angular_distance` at :math:`t = 0`.
"""

initial_angular_speed = clone_as_symbol(symbols.angular_speed, subscript="0")
"""
:symbols:`angular_speed` at :math:`t = 0`.
"""

angular_acceleration = symbols.angular_acceleration
"""
Constant :symbols:`angular_acceleration`.
"""

time = symbols.time
"""
:symbols:`time` at which :attr:`~final_angular_position` is measured.
"""

law = Eq(
    final_angular_position,
    initial_angular_position + initial_angular_speed * time + angular_acceleration * time**2 / 2,
)
"""
:laws:symbol::

:laws:latex::
"""

# Derive law from definitions of angular velocity and acceleration

_angular_velocity_formula = dsolve(
    angular_acceleration_def.definition.subs(angular_acceleration_def.time, time),
    angular_acceleration_def.angular_speed(time),
).rhs.subs(
    angular_acceleration_def.angular_acceleration(time),
    angular_acceleration,
).doit()

_angular_velocity = symbols.angular_speed
_angular_velocity_derived = solve(
    [
    Eq(initial_angular_speed, _angular_velocity_formula.subs(time, 0)),
    Eq(_angular_velocity, _angular_velocity_formula)
    ],
    ("C1", _angular_velocity),
    dict=True,
)[0][_angular_velocity]

_angular_displacement_formula = dsolve(
    angular_velocity_def.definition.subs(angular_velocity_def.time, time),
    angular_velocity_def.angular_distance(time),
).rhs.subs(
    angular_velocity_def.angular_speed(time),
    _angular_velocity_derived,
).doit()

_angular_displacement_derived = solve(
    [
    # initial angular displacement is 0 by condition
    Eq(initial_angular_position, _angular_displacement_formula.subs(time, 0)),
    Eq(final_angular_position, _angular_displacement_formula)
    ],
    ("C1", final_angular_position),
    dict=True,
)[0][final_angular_position]

assert expr_equals(_angular_displacement_derived, law.rhs)


@validate_input(
    initial_angular_position_=initial_angular_position,
    initial_angular_velocity_=initial_angular_speed,
    angular_acceleration_=angular_acceleration,
    time_=time,
)
@validate_output(final_angular_position)
def calculate_angular_displacement(
    initial_angular_position_: Quantity | float,
    initial_angular_velocity_: Quantity,
    angular_acceleration_: Quantity,
    time_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        initial_angular_position: scale_factor(initial_angular_position_),
        initial_angular_speed: initial_angular_velocity_,
        angular_acceleration: angular_acceleration_,
        time: time_,
    })
    return Quantity(result)

r"""
Angular position via constant angular speed and time
====================================================

When a body is rotating around a fixed axis with a constant angular speed, its angular
position is a linear function of time.

#. The axis is fixed.
#. Angular speed is constant, i.e. :math:`\frac{d \omega}{d t} = 0.`
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
from symplyphysics.definitions import angular_speed_is_angular_distance_derivative as angular_velocity_def

final_angular_position = symbols.angular_distance
"""
:symbols:`angular_distance` at :attr:`~time`.
"""

initial_angular_position = clone_as_symbol(symbols.angular_distance, subscript="0")
"""
:symbols:`angular_distance` at :math:`t = 0`.
"""

angular_speed = symbols.angular_speed
"""
Constant :symbols:`angular_speed`.
"""

time = symbols.time
"""
:symbols:`time` at which :attr:`~final_angular_position` is measured.
"""

law = Eq(
    final_angular_position,
    initial_angular_position + angular_speed * time,
)
"""
:laws:symbol::

:laws:latex::
"""

# Derive law from definition of angular velocity

_angular_position_formula = dsolve(
    angular_velocity_def.definition.subs(angular_velocity_def.time, time),
    angular_velocity_def.angular_distance(time),
).rhs.subs(
    angular_velocity_def.angular_speed(time),
    angular_speed,
).doit()

_c1 = solve(Eq(initial_angular_position, _angular_position_formula.subs(time, 0)), "C1")[0]

_angular_position_derived = _angular_position_formula.subs("C1", _c1)

assert expr_equals(_angular_position_derived, law.rhs)


@validate_input(
    initial_angular_position_=initial_angular_position,
    angular_velocity_=angular_speed,
    time_=time,
)
@validate_output(final_angular_position)
def calculate_angular_position(
    initial_angular_position_: Quantity | float,
    angular_velocity_: Quantity,
    time_: Quantity,
) -> Quantity:
    initial_angular_position_ = (initial_angular_position_.scale_factor if isinstance(
        initial_angular_position_, Quantity) else initial_angular_position_)
    result = law.rhs.subs({
        initial_angular_position: initial_angular_position_,
        angular_speed: angular_velocity_,
        time: time_,
    })
    return Quantity(result)

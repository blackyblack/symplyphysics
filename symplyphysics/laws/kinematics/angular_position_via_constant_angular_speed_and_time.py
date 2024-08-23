r"""
Angular position via constant angular speed and time
====================================================

When a body is rotating around a fixed axis with a constant angular speed, its angular
position is a linear function of time.

#. The axis is fixed.
#. :math:`\omega = \text{const} \Leftrightarrow \frac{d \omega}{d t} = 0.`
"""

from sympy import Eq, solve, dsolve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    angle_type,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import angular_speed_is_angular_distance_derivative as angular_velocity_def

final_angular_position = Symbol("final_angular_position", angle_type)
r"""
Angular position at time :math:`t`.

Symbol:
    :code:`theta`

Latex:
    :math:`\theta`
"""

initial_angular_position = Symbol("initial_angular_position", angle_type)
r"""
Angular position at :math:`t = 0`.

Symbol:
    :code:`theta_0`

Latex:
    :math:`\theta_0`
"""

angular_speed = Symbol("angular_speed", angle_type / units.time)
r"""
Constant angular speed.

Symbol:
    :code:`w`

Latex:
    :math:`\omega`
"""

time = Symbol("time", units.time)
r"""
Time at which :math:`\theta` is measured.

Symbol:
    :code:`t`
"""

law = Eq(
    final_angular_position,
    initial_angular_position + angular_speed * time,
)
r"""
:code:`theta = theta_0 + w * t`

Latex:
    .. math::
        \theta = \theta_0 + \omega t
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

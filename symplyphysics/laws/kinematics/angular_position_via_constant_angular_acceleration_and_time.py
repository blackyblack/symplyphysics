r"""
Angular position via constant angular acceleration and time
===========================================================

If a body is rotating with a constant acceleration, its angular position is a quadratic function of time.

#. The axis is fixed.
#. Angular acceleration is constant, i.e. :math:`\frac{d \alpha}{d t} = 0.`
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
from symplyphysics.core.symbols.quantities import scale_factor
from symplyphysics.definitions import (
    angular_acceleration_is_angular_speed_derivative as angular_acceleration_def,
    angular_speed_is_angular_distance_derivative as angular_velocity_def,
)

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

initial_angular_speed = Symbol("initial_angular_speed", angle_type / units.time)
r"""
Angular speed at :math:`t = 0`.

Symbol:
    :code:`w_0`

Latex:
    :math:`\omega_0`
"""

angular_acceleration = Symbol("angular_acceleration", angle_type / units.time**2)
r"""
Constant angular acceleration.

Symbol:
    :code:`alpha`

Latex:
    :math:`\alpha`
"""

time = Symbol("time", units.time)
r"""
Time at which :math:`\theta` is measured.

Symbol:
    :code:`t`
"""

law = Eq(
    final_angular_position,
    initial_angular_position + initial_angular_speed * time + angular_acceleration * time**2 / 2,
)
r"""
:code:`theta = theta_0 + w_0 * t + 1/2 * alpha * t^2`

Latex:
    .. math::
        \theta = \theta_0 + \omega_0 t + \frac{1}{2} \alpha t^2
"""


# Derive law from definitions of angular velocity and acceleration

_angular_velocity_formula = dsolve(
    angular_acceleration_def.definition.subs(angular_acceleration_def.time, time),
    angular_acceleration_def.angular_speed(time),
).rhs.subs(
    angular_acceleration_def.angular_acceleration(time),
    angular_acceleration,
).doit()

_angular_velocity = Symbol("_angular_velocity", angle_type / units.time)
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

r"""
Angular speed via constant angular acceleration and time
========================================================

If a body is rotating around a fixed axis with constant angular acceleration, its angular
speed is a linear function of time.

**Conditions:**

#. The axis is fixed.
#. Angular acceleration is constant, i.e. :math:`\frac{d \alpha}{d t} = 0.`
"""

from sympy import Eq, dsolve, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    angle_type,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import (
    angular_acceleration_is_angular_speed_derivative as angular_acceleration_def,)

final_angular_speed = Symbol("final_angular_speed", angle_type / units.time)
r"""
Angular speed at time :math:`t`.

Symbol:
    :code:`w`

Latex:
    :math:`\omega`
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
Time at which :math:`\omega` is measured.

Symbol:
    :code:`t`
"""

law = Eq(final_angular_speed, initial_angular_speed + angular_acceleration * time)
r"""
:code:`w = w_0 + alpha * t`

Latex:
    .. math::
        \omega = \omega_0 + \alpha t
"""

# Derive this law from definition of angular acceleration

_angular_velocity_formula = dsolve(
    angular_acceleration_def.definition.subs(angular_acceleration_def.time, time),
    angular_acceleration_def.angular_speed(time),
).rhs.subs(
    angular_acceleration_def.angular_acceleration(time),
    angular_acceleration,
).doit()

_angular_velocity_derived = solve([
    Eq(initial_angular_speed, _angular_velocity_formula.subs(time, 0)),
    Eq(final_angular_speed, _angular_velocity_formula)
], ("C1", final_angular_speed),
    dict=True)[0][final_angular_speed]

assert expr_equals(_angular_velocity_derived, law.rhs)


@validate_input(
    initial_angular_velocity_=initial_angular_speed,
    angular_acceleration_=angular_acceleration,
    time_=time,
)
@validate_output(final_angular_speed)
def calculate_angular_velocity(
    initial_angular_velocity_: Quantity,
    angular_acceleration_: Quantity,
    time_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        initial_angular_speed: initial_angular_velocity_,
        angular_acceleration: angular_acceleration_,
        time: time_,
    })
    return Quantity(result)

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
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import (
    angular_acceleration_is_angular_speed_derivative as angular_acceleration_def,)

final_angular_speed = symbols.angular_speed
"""
:symbols:`angular_speed` at :attr:`~time`.
"""

initial_angular_speed = clone_as_symbol(symbols.angular_speed, subscript="0")
"""
:symbols:`angular_speed` at :math:`t = 0`.
"""

angular_acceleration = symbols.angular_acceleration
r"""
Constant :symbols:`angular_acceleration`.
"""

time = symbols.time
"""
:symbols:`time` at which :attr:`~final_angular_speed` is measured.
"""

law = Eq(final_angular_speed, initial_angular_speed + angular_acceleration * time)
"""
:laws:symbol::

:laws:latex::
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

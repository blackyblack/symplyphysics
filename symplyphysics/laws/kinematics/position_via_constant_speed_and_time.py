r"""
Position via constant speed and time
====================================

If a body is moving with constant speed, its position can be expressed as a linear function
of time.

**Conditions:**

#. :math:`v = \text{const} \Leftrightarrow \frac{d v}{d t} = 0.`
"""

from sympy import (Eq, solve, dsolve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import speed_is_distance_derivative as velocity_definition

final_position = Symbol("distance_function", units.length)
r"""
Position at time :math:`t`.

Symbol:
    :code:`x`
"""

initial_position = Symbol("initial_position", units.length)
r"""
Position at :math:`t = 0`.

Symbol:
    :code:`x_0`

Latex:
    :math:`x_0`
"""

speed = Symbol("speed", units.velocity, constant=True)
"""
Constant speed.

Symbol:
    :code:`v`
"""

time = Symbol("time", units.time)
r"""
Time at which :math:`x` is measured.

Symbol:
    :code:`t`
"""

law = Eq(final_position, initial_position + speed * time)
r"""
:code:`x = x_0 + v * t`

Latex:
    .. math::
        x = x_0 + v t
"""

# Derive the same law from velocity definition

_constant_velocity_movement_definition = velocity_definition.definition.subs({
    velocity_definition.speed(velocity_definition.time): speed,
    velocity_definition.time: time
})
_dsolved_movement = dsolve(_constant_velocity_movement_definition,
    velocity_definition.distance(time))

# Prove that derived movement function equals to law.rhs, given C1 = initial_position
assert (expr_equals(_dsolved_movement.rhs.subs("C1", initial_position), law.rhs))


@validate_input(initial_distance_=initial_position,
    velocity_=speed,
    time_=time)
@validate_output(final_position)
def calculate_distance(initial_distance_: Quantity, velocity_: Quantity,
    time_: Quantity) -> Quantity:
    result_expr = solve(law, final_position, dict=True)[0][final_position]
    result_expr_substituted = result_expr.subs({
        initial_position: initial_distance_,
        speed: velocity_,
        time: time_
    })
    return Quantity(result_expr_substituted)

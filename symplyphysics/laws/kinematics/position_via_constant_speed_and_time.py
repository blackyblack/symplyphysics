r"""
Position via constant speed and time
====================================

If a body is moving with constant speed, its position can be expressed as a linear function
of time.

**Conditions:**

#. Speed is constant, i.e. :math:`\frac{d v}{d t} = 0.`

**Links:**

#. `Wikipedia, derivable from the vector counterpart of this law <https://en.wikipedia.org/wiki/Kinematics#Relative_acceleration>`__.
"""

from sympy import (Eq, solve, dsolve)
from symplyphysics import (Quantity, validate_input, validate_output, symbols, clone_as_symbol)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import speed_is_distance_derivative as velocity_definition

final_position = symbols.position
"""
:symbols:`position` at :attr:`~time`.
"""

initial_position = clone_as_symbol(symbols.position, subscript="0")
"""
:symbols:`position` at :math:`t = 0`.
"""

speed = symbols.speed
"""
Constant :symbols:`speed`.
"""

time = symbols.time
"""
:symbols:`time` at which :attr:`~final_position` is measured.
"""

law = Eq(final_position, initial_position + speed * time)
"""
:laws:symbol::

:laws:latex::
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


@validate_input(initial_distance_=initial_position, velocity_=speed, time_=time)
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

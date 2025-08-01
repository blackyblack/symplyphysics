"""
Position via constant acceleration and time
===========================================

If a body is moving with a constant acceleration, its position in space is a quadratic function of time.

**Conditions:**

#. Acceleration is constant, i.e. :math:`\\frac{d a}{d t} = 0.`

**Links:**

#. `Wikipedia, vector counterpart of this law <https://en.wikipedia.org/wiki/Kinematics#Relative_acceleration>`__.

..
    TODO: make a vector counterpart of this law
"""

from sympy import Eq, solve, dsolve
from symplyphysics import symbols, Quantity, validate_input, validate_output, clone_as_symbol
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import speed_is_distance_derivative as _velocity_definition
from symplyphysics.definitions import acceleration_is_speed_derivative as _acceleration_definition

final_position = symbols.position
"""
:symbols:`position` at :attr:`~time`.
"""

initial_position = clone_as_symbol(symbols.position, subscript="0")
"""
:symbols:`position` at :math:`t = 0`.
"""

initial_speed = clone_as_symbol(symbols.speed, subscript="0")
"""
:symbols:`speed` at :math:`t = 0`. Note that it is the *projection* of the velocity vector on the
axis of the body's movement.
"""

acceleration = symbols.acceleration
"""
Constant :symbols:`acceleration`. Note that it is the *projection* of the acceleration vector on the
axis of the body's movement.
"""

time = symbols.time
"""
:symbols:`time` at which :attr:`~final_position` is measured.
"""

law = Eq(final_position, initial_position + initial_speed * time + acceleration * time**2 / 2)
"""
:laws:symbol::

:laws:latex::
"""

# Derive the same law from velocity and acceleration definitions

_constant_acceleration_definition = _acceleration_definition.definition.subs({
    _acceleration_definition.acceleration(_acceleration_definition.time): acceleration,
    _acceleration_definition.time: time
})
_dsolved_velocity = dsolve(_constant_acceleration_definition, _acceleration_definition.speed(time))
_constant_accelerated_velocity_function = _dsolved_velocity.rhs

_constant_accelerated_movement_definition = _velocity_definition.definition.subs({
    _velocity_definition.speed(_velocity_definition.time): _constant_accelerated_velocity_function,
    _velocity_definition.time: time
})
_dsolved_movement = dsolve(_constant_accelerated_movement_definition,
    _velocity_definition.distance(time))
_constant_accelerated_movement_function = _dsolved_movement.rhs

_derived_law = Eq(initial_position, _constant_accelerated_movement_function)

# Prove that _constant_accelerated_movement_function equals to law.rhs, given C1 = initial_speed,
# C2 = initial initial_position = 0
assert expr_equals(_derived_law.rhs.subs({"C1": initial_speed, "C2": initial_position}), law.rhs)


@validate_input(initial_position_=initial_position,
    initial_velocity_=initial_speed,
    acceleration_=acceleration,
    time_=time)
@validate_output(final_position)
def calculate_distance(initial_position_: Quantity, initial_velocity_: Quantity,
    acceleration_: Quantity, time_: Quantity) -> Quantity:
    result_expr = solve(law, final_position, dict=True)[0][final_position]
    result_expr_substituted = result_expr.subs({
        initial_position: initial_position_,
        initial_speed: initial_velocity_,
        acceleration: acceleration_,
        time: time_
    })
    return Quantity(result_expr_substituted)

"""
Maximum height from initial speed
=================================

The maximum height that a body thrown vertically will rise to depends on the initial speed.

**Notes:**

#. If the body is launched at an arbitrary angle to the horizontal, this formula is still
   applicable but you must use the projection of the velocity vector on the :math:`z` axis as the
   initial speed.

**Conditions:**

#. The body is thrown "upwards", i.e. in the direction opposite to that of the free fall
   acceleration vector.

**Links:**

#. `Wikipedia, fifth formula <https://en.wikipedia.org/wiki/Projectile_motion#Maximum_height_of_projectile>`__.

..
    TODO: move to kinetimatics
"""

from sympy import Eq, solve, pi
from symplyphysics import (Quantity, validate_input, validate_output, symbols, quantities,
    clone_as_symbol)
from symplyphysics.core.expr_comparisons import expr_equals

from symplyphysics.laws.geometry import (
    scalar_projection_is_vector_length_times_cosine_of_angle as _projection_law,)
from symplyphysics.laws.kinematics import (
    speed_via_constant_acceleration_and_time as _speed_law,
    position_via_constant_acceleration_and_time as _position_law,
)

maximum_height = symbols.height
"""
The maximum :symbols:`height` that the object will reach.
"""

initial_speed = clone_as_symbol(symbols.speed)
"""
The initial :symbols:`speed` of the object. Note that it is the *projection* of the body's velocity
vector on the :math:`z`-axis.
"""

law = Eq(maximum_height, initial_speed**2 / (2 * quantities.acceleration_due_to_gravity))
"""
:laws:symbol::

:laws:latex::
"""

### Derive law

_acceleration = _projection_law.law.rhs.subs({
    _projection_law.vector_length: quantities.acceleration_due_to_gravity,
    _projection_law.angle: pi
})

_speed_eqn = _speed_law.law.subs({
    _speed_law.final_speed: 0,  # the body exhausts all its kinetic energy and comes to a stop
    _speed_law.initial_speed: initial_speed,
    _speed_law.acceleration: _acceleration,
})

_time = _speed_law.time

_position_eqn = _position_law.law.subs({
    _position_law.final_position: maximum_height,
    _position_law.initial_position:
        0,  # `maximum_height` is the relative distance along the z-axis the body travels during the flight
    _position_law.initial_speed: initial_speed,
    _position_law.acceleration: _acceleration,
    _position_law.time: _time,
})

_maximum_height_derived = solve(
    (_speed_eqn, _position_eqn),
    (maximum_height, _time),
    dict=True,
)[0][maximum_height]

assert expr_equals(_maximum_height_derived, law.rhs)

### Let's show that body must be thrown upwards, i.e. the initial speed projection must be positive

# Suppose the initial speed projection is negative:
_initial_speed_neg = clone_as_symbol(initial_speed, negative=True)

# Time is positive: the body will attain the maximum height at a later time than when it is thrown
_time_pos = clone_as_symbol(_time, positive=True)

_speed_eqn_subs = _speed_eqn.subs({
    initial_speed: _initial_speed_neg,
    _time: _time_pos,
})

# The equation above isn't solvable
assert not solve(_speed_eqn_subs, _time_pos)


@validate_input(initial_velocity_=initial_speed)
@validate_output(maximum_height)
def calculate_maximum_height(initial_velocity_: Quantity) -> Quantity:
    result_maximum_height = solve(law, maximum_height, dict=True)[0][maximum_height]
    result_expr = result_maximum_height.subs({initial_speed: initial_velocity_})
    return Quantity(result_expr)


# UNIQUE_LAW_ID: 195

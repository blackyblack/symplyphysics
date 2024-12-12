r"""
Tangential acceleration via angular acceleration and radius
===========================================================

The tangential acceleration of a rotating body represents the change in magnitude
of the velocity vector, and its vector is tangent to the path of the body.

**Conditions:**

#. Radius is constant, i.e. :math:`\frac{d r}{d t} = 0.`

**Links:**

#. Equation 10-22 on p. 269 of "Fundamentals of Physics" by David Halladay et al., 10th Ed.
"""

from sympy import Eq, solve, Derivative
from symplyphysics import (
    clone_as_symbol,
    symbols,
    Quantity,
    validate_input,
    validate_output,
    clone_as_function,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.kinematics import speed_via_angular_speed_and_radius as linear_velocity_law
from symplyphysics.definitions import angular_acceleration_is_angular_speed_derivative as angular_acceleration_def
from symplyphysics.definitions import acceleration_is_speed_derivative as acceleration_def

tangential_acceleration = clone_as_symbol(symbols.acceleration,
    display_symbol="a_t",
    display_latex="a_\\tau")
"""
Tangential :symbols:`acceleration`.
"""

angular_acceleration = symbols.angular_acceleration
"""
:symbols:`angular_acceleration`.
"""

radius_of_curvature = symbols.radius_of_curvature
"""
Instantaneous :symbols:`radius_of_curvature`.
"""

law = Eq(tangential_acceleration, angular_acceleration * radius_of_curvature)
"""
:laws:symbol::

:laws:latex::
"""

# Derive the law from the law of linear velocity in case of a rotating body
# Condition: radius of rotation remains constant.

_time = symbols.time
_linear_velocity = clone_as_function(symbols.speed, [_time])
_angular_velocity = clone_as_function(symbols.angular_speed, [_time])

_linear_velocity_law_sub = linear_velocity_law.law.subs({
    linear_velocity_law.speed: _linear_velocity(_time),
    linear_velocity_law.angular_speed: _angular_velocity(_time),
    linear_velocity_law.radius_of_curvature: radius_of_curvature,
})

# Differentiate both sides of the equation w.r.t. _time
_diff_linear_velocity_law = Eq(_linear_velocity_law_sub.lhs.diff(_time),
    _linear_velocity_law_sub.rhs.diff(_time))

# alpha = d(omega)/dt
_angular_acceleration_def_sub = angular_acceleration_def.definition.subs(
    angular_acceleration_def.time, _time)
_angular_acceleration_def_sub = _angular_acceleration_def_sub.subs(
    angular_acceleration_def.angular_speed(_time), _angular_velocity(_time))

_linear_velocity_derivative = solve([_diff_linear_velocity_law, _angular_acceleration_def_sub],
    (Derivative(_angular_velocity(_time), _time), Derivative(_linear_velocity(_time), _time)),
    dict=True)[0][Derivative(_linear_velocity(_time), _time)]
linear_velocity_derivative_eq = Eq(Derivative(_linear_velocity(_time), _time),
    _linear_velocity_derivative)

# a_t = dv/dt
_acceleration_def_sub = acceleration_def.definition.subs(acceleration_def.time, _time)
_acceleration_def_sub = _acceleration_def_sub.subs(acceleration_def.speed(_time),
    _linear_velocity(_time))
_tangential_acceleration_value = solve([linear_velocity_derivative_eq, _acceleration_def_sub],
    (Derivative(_linear_velocity(_time), _time), acceleration_def.acceleration(_time)),
    dict=True)[0][acceleration_def.acceleration(_time)]
_tangential_acceleration_derived = _tangential_acceleration_value.subs(
    angular_acceleration_def.angular_acceleration(_time), angular_acceleration)

assert expr_equals(
    _tangential_acceleration_derived,
    law.rhs,
)


@validate_input(
    angular_acceleration_=angular_acceleration,
    rotation_radius_=radius_of_curvature,
)
@validate_output(tangential_acceleration)
def calculate_tangential_acceleration(angular_acceleration_: Quantity,
    rotation_radius_: Quantity) -> Quantity:
    result_expr = solve(law, tangential_acceleration)[0]
    result = result_expr.subs({
        angular_acceleration: angular_acceleration_,
        radius_of_curvature: rotation_radius_,
    })
    return Quantity(result)

r"""
Tangential acceleration via angular acceleration and radius
===========================================================

The tangential acceleration of a rotating body represents the change in magnitude
of the velocity vector, and its vector is tangent to the path of the body.

**Conditions:**

#. :math:`r = \text{const} \Leftrightarrow = \frac{d r}{d t} = 0.`
"""

from sympy import Eq, solve, Derivative
from symplyphysics import (clone_symbol, symbols, units, Quantity, Symbol, Function
    , validate_input, validate_output, angle_type)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.kinematics import speed_via_angular_speed_and_radius as linear_velocity_law
from symplyphysics.definitions import angular_acceleration_is_angular_speed_derivative as angular_acceleration_def
from symplyphysics.definitions import acceleration_is_speed_derivative as acceleration_def

tangential_acceleration = clone_symbol(symbols.kinematics.acceleration, "tangential_acceleration")
r"""
Tangential :attr:`~symplyphysics.symbols.kinematics.acceleration`.

Symbol:
    :code:`a_t`

Latex:
    :math:`a_\tau`
"""

angular_acceleration = Symbol("angular_acceleration", angle_type / units.time**2)
r"""
Angular acceleration.

Symbol:
    :code:`alpha`

Latex:
    :math:`\alpha`
"""

radius_of_curvature = Symbol("radius_of_curvature", units.length)
"""
Instantaneous radius of curvature.

Symbol:
    :code:`r`
"""

law = Eq(tangential_acceleration, angular_acceleration * radius_of_curvature)
r"""
:code:`a_t = alpha * r`

Latex:
    .. math::
        a_\tau = \alpha r
"""

# Derive the law from the law of linear velocity in case of a rotating body
# Condition: radius of rotation remains constant.

_linear_velocity = Function("_linear_velocity", units.velocity)
_angular_velocity = Function("_angular_velocity", angle_type / units.time)
_time = Symbol("_time", units.time)

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

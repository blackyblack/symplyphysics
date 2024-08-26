"""
Angular position is arc length over radius
==========================================

To describe the rotation of a rigid body about a fixed rotational axis, a reference line is
assumed to be fixed in the body, perpendicular to that axis and rotating with the body. The angular
position of this is measured relative to a fixed direction and is expressed as a ratio of the arc
length of a circular path and its radius (distance to the axis).
"""

from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    angle_type,
)

angular_position = Symbol("angular_position", angle_type)
r"""
Angular position of the body.

Symbol:
    :code:`theta`

Latex:
    :math:`\theta`
"""

arc_length = Symbol("arc_length", units.length)
"""
Arc length.

Symbol:
    :code:`s`
"""

distance_to_axis = Symbol("distance_to_axis", units.length)
"""
Distance to rotational axis, or radius of rotation.

Symbol:
    :code:`r`
"""

law = Eq(angular_position, arc_length / distance_to_axis)
r"""
:code:`theta = s / r`

Latex:
    .. math::
        \theta = \frac{s}{r}
"""


@validate_input(arc_length_=arc_length, path_radius_=distance_to_axis)
@validate_output(angular_position)
def calculate_angular_position(arc_length_: Quantity, path_radius_: Quantity) -> Quantity:
    result = law.rhs.subs({
        arc_length: arc_length_,
        distance_to_axis: path_radius_,
    })
    result_angle = Quantity(result)
    return result_angle

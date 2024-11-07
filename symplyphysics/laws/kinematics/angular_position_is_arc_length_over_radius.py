"""
Angular position is arc length over radius
==========================================

To describe the rotation of a rigid body about a fixed rotational axis, a reference line is
assumed to be fixed in the body, perpendicular to that axis and rotating with the body. The angular
position of this is measured relative to a fixed direction and is expressed as a ratio of the arc
length of a circular path and its radius (distance to the axis).
"""

from sympy import Eq
from symplyphysics import Quantity, validate_input, validate_output, symbols

angular_position = symbols.angular_distance
r"""
:symbols:`angular_distance` of the body.
"""

arc_length = symbols.arc_length
"""
:symbols:`arc_length` of the curve traced by the body's movement.
"""

distance_to_axis = symbols.distance_to_axis
"""
:symbols:`distance_to_axis` of rotation, or radius of rotation.
"""

law = Eq(angular_position, arc_length / distance_to_axis)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(arc_length_=arc_length, path_radius_=distance_to_axis)
@validate_output(angular_position)
def calculate_angular_position(arc_length_: Quantity, path_radius_: Quantity) -> Quantity:
    result = law.rhs.subs({
        arc_length: arc_length_,
        distance_to_axis: path_radius_,
    })
    return Quantity(result)

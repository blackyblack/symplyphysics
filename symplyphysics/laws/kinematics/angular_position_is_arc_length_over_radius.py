"""
Angular position is arc length over radius
==========================================

To describe the rotation of a rigid body about a fixed rotational axis, a reference line is
assumed to be fixed in the body, perpendicular to that axis and rotating with the body. The angular
position of this is measured relative to a fixed direction and is expressed as a ratio of the arc
length of a circular path and its radius (distance to the axis).

**Links:**

#. `Wikipedia, similar concept <https://en.wikipedia.org/wiki/Central_angle#Formulas>`__.

#. `openstax, table 6.2, first line <https://openstax.org/books/physics/pages/6-3-rotational-motion>`__.
"""

from sympy import Eq, solve
from symplyphysics import Quantity, validate_input, validate_output, symbols, clone_as_symbol
from symplyphysics.core.expr_comparisons import expr_equals

from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, CoordinateVector
from symplyphysics.core.experimental.vectors import VectorNorm
from symplyphysics.core.experimental.solvers import solve_for_vector

from symplyphysics.laws.kinematics.vector import (
    displacement_is_angular_displacement_cross_radius as _linear_displacement_law,)

angular_position = clone_as_symbol(symbols.angular_distance, positive=True)
"""
:symbols:`angular_distance` of the body.
"""

arc_length = clone_as_symbol(symbols.arc_length, positive=True)
"""
:symbols:`arc_length` of the curve traced by the body's movement.
"""

distance_to_axis = clone_as_symbol(symbols.distance_to_axis, positive=True)
"""
:symbols:`distance_to_axis` of rotation, or radius of rotation.
"""

law = Eq(angular_position, arc_length / distance_to_axis)
"""
:laws:symbol::

:laws:latex::
"""

# Derive from vector law

_angular_displacement = CoordinateVector([angular_position, 0, 0], CARTESIAN)
_radius_vector = CoordinateVector([0, distance_to_axis, 0], CARTESIAN)

_linear_displacement_expr = solve_for_vector(
    _linear_displacement_law.law,
    _linear_displacement_law.linear_displacement,
).subs({
    _linear_displacement_law.angular_displacement: _angular_displacement,
    _linear_displacement_law.rotation_radius_vector: _radius_vector,
}).doit()

_arc_length_derived = VectorNorm(_linear_displacement_expr).doit()

_arc_length_expected = solve(law, arc_length)[0]

assert expr_equals(_arc_length_derived, _arc_length_expected)


@validate_input(arc_length_=arc_length, path_radius_=distance_to_axis)
@validate_output(angular_position)
def calculate_angular_position(arc_length_: Quantity, path_radius_: Quantity) -> Quantity:
    result = law.rhs.subs({
        arc_length: arc_length_,
        distance_to_axis: path_radius_,
    })
    return Quantity(result)


# UNIQUE_LAW_ID: 433

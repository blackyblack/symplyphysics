"""
Linear displacement is angular displacement cross radius
========================================================

Assuming a body rotating around a fixed axis, the vector of its linear displacement can be expressed
as the cross product of the pseudovector of angular displacement and the radius vector of rotation.

**Conditions:**

#. The axis is fixed.
#. Angular displacement pseudovector and radius vector must be orthogonal to one another.

**Links:**

#. `Physics LibreTexts, formula 11.1.4 <https://phys.libretexts.org/Bookshelves/University_Physics/Book%3A_Introductory_Physics_-_Building_Models_to_Describe_Our_World_(Martin_Neary_Rinaldo_and_Woodman)/11%3A_Rotational_dynamics/11.01%3A_Rotational_kinematic_vectors>`__.
"""

from sympy import Eq
from symplyphysics import Quantity, validate_input, validate_output, symbols

from symplyphysics.core.approx import approx_equal_numbers
from symplyphysics.core.experimental.vectors import clone_as_vector_symbol, VectorCross, VectorDot
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

linear_displacement = clone_as_vector_symbol(symbols.distance)
"""
Vector of the body's linear displacement. See :symbols:`distance`.
"""

angular_displacement = clone_as_vector_symbol(symbols.angular_distance)
"""
Pseudovector of the body's angular displacement. See :symbols:`angular_distance`. It is parallel
to the rotation axis.
"""

rotation_radius_vector = clone_as_vector_symbol(symbols.distance_to_axis)
"""
Radius vector pointing away from the rotation axis perpendicular to it. See
:symbols:`distance_to_axis`.
"""

law = Eq(
    linear_displacement,
    VectorCross(angular_displacement, rotation_radius_vector),
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    angular_displacement_=angular_displacement,
    rotation_radius_=rotation_radius_vector,
)
@validate_output(linear_displacement)
def calculate_linear_displacement(
    angular_displacement_: QuantityCoordinateVector,
    rotation_radius_: QuantityCoordinateVector,
) -> QuantityCoordinateVector:
    dot = Quantity(VectorDot(angular_displacement_, rotation_radius_))

    if not approx_equal_numbers(dot.scale_factor, 0, absolute_tolerance=1e-10):
        raise ValueError("Angular displacement and rotation radius vector must be orthogonal")

    result = law.rhs.subs({
        angular_displacement: angular_displacement_,
        rotation_radius_vector: rotation_radius_,
    })

    return QuantityCoordinateVector.from_expr(result)


# UNIQUE_LAW_ID: 453

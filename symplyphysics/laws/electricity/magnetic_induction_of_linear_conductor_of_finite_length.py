"""
Magnetic flux density of linear conductor of finite length
==========================================================

Let there be a rectilinear conductor of finite length. Then its magnetic flux density
will depend on the magnitude of the current and the material. It also depends on the
perpendicular distance to the conductor and on the angles between the lines drawn from
the ends of the conductor to the point and the conductor.

**Conditions:**

#. Conductor should be rectilinear.
#. Length of the conductor is finite.

**Links:**

#. `Physics LibreTexts — infinite wire case <https://phys.libretexts.org/Bookshelves/University_Physics/University_Physics_(OpenStax)/University_Physics_II_-_Thermodynamics_Electricity_and_Magnetism_(OpenStax)/12%3A_Sources_of_Magnetic_Fields/12.03%3A_Magnetic_Field_due_to_a_Thin_Straight_Wire>`__.
"""

from sympy import Eq, solve, pi, cos, simplify, integrate, cot, refine, Q, sin
from symplyphysics import (Quantity, validate_input, validate_output, symbols, clone_as_symbol)
from symplyphysics.laws.electricity.vector import magnetic_field_due_to_constant_filamentary_current as _biot_savart_law

from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.experimental.vectors import VectorNorm, VectorCross, VectorDot
from symplyphysics.core.experimental.solvers import vector_equals
from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, CoordinateVector

magnetic_flux_density = symbols.magnetic_flux_density
"""
:symbols:`magnetic_flux_density` through the conductor.
"""

absolute_permeability = symbols.absolute_permeability
"""
:symbols:`absolute_permeability` of the medium.
"""

current = symbols.current
"""
:symbols:`current` running through the conductor.
"""

first_angle = clone_as_symbol(symbols.angle, subscript="1")
"""
:symbols:`angle` between origin and the first end of the conductor.
"""

second_angle = clone_as_symbol(symbols.angle, subscript="2")
"""
:symbols:`angle` between origin and the second end of the conductor.
"""

radial_distance = clone_as_symbol(symbols.euclidean_distance, positive=True)
"""
Perpendicular :symbols:`euclidean_distance` to the conductor.
"""

law = Eq(
    magnetic_flux_density,
    absolute_permeability * current * (cos(first_angle) + cos(second_angle)) /
    (4 * pi * radial_distance),
)
"""
:laws:symbol::

:laws:latex::
"""

# Derive from the Biot-Savart law
# Link: <https://phys.libretexts.org/Bookshelves/University_Physics/University_Physics_(OpenStax)/University_Physics_II_-_Thermodynamics_Electricity_and_Magnetism_(OpenStax)/12%3A_Sources_of_Magnetic_Fields/12.03%3A_Magnetic_Field_due_to_a_Thin_Straight_Wire>

# Let `P` be the point where the magnetic flux density is measured, `PO` be the line segment
# perpendicular to the wire, i.e. `|PO|` is the radial distance to the wire. Let the `x`-axis run
# in the (conventional) direction of the current and the `y`-axis, which contains `P`, is
# orthogonal to it. Then the infinitesimal displacement of the contour element also lies on the
# `x`-axis.

_position_vector = CoordinateVector([0, radial_distance, 0], CARTESIAN)

_contour_element_position = clone_as_symbol(symbols.position, real=True)
_contour_element_length = clone_as_symbol(
    symbols.euclidean_distance,
    display_symbol="dx",
    positive=True,
)

_contour_element_position_vector = CoordinateVector([_contour_element_position, 0, 0], CARTESIAN)
_contour_element_displacement = CoordinateVector([_contour_element_length, 0, 0], CARTESIAN)

_magnetic_flux_density_change_vector = _biot_savart_law.law.rhs.subs({
    _biot_savart_law.absolute_permeability: absolute_permeability,
    _biot_savart_law.current: current,
    _biot_savart_law.position_vector: _position_vector,
    _biot_savart_law.contour_element_position_vector: _contour_element_position_vector,
    _biot_savart_law.contour_element_displacement: _contour_element_displacement
}).doit()

# Convert linear combinations of vectors in a single vector
for norm in _magnetic_flux_density_change_vector.atoms(VectorNorm):
    _magnetic_flux_density_change_vector = _magnetic_flux_density_change_vector.subs(
        norm, VectorNorm(CoordinateVector.from_expr(*norm.args)))
_magnetic_flux_density_change_vector = CoordinateVector.from_expr(
    _magnetic_flux_density_change_vector)

# Check that the vector of magnetic flux density change lies parallel to the `z`-axis
_e_z = CoordinateVector([0, 0, 1], CARTESIAN)
assert vector_equals(VectorCross(_magnetic_flux_density_change_vector, _e_z), 0)

_magnetic_flux_density_change_z_component = simplify(
    VectorDot(_magnetic_flux_density_change_vector, _e_z))

# This is derivable from the geometry of the system (see link for the figure). Note that angles are
# measured as the inner angles, but after integration we get angles that are measured between the
# position vector from the contour element and the negative direction of the `x`-axis. That is, if
# the first angle is acute, then the first ("left") end of the wire is located to the "left" (in
# the negative half-ray) of the `x`-axis, therefore we have to take the complement of that angle to
# use in the integration formula.
_first_position = radial_distance * cot(pi - first_angle)
_second_position = radial_distance * cot(second_angle)

_total_magnetic_flux_density = integrate(
    _magnetic_flux_density_change_z_component.subs(_contour_element_length, 1),
    (_contour_element_position, _first_position, _second_position),
)
_total_magnetic_flux_density = refine(
    simplify(_total_magnetic_flux_density),
    Q.positive(sin(first_angle)) & Q.positive(sin(second_angle)),
)
_total_magnetic_flux_density = simplify(_total_magnetic_flux_density)

assert expr_equals(law.rhs, _total_magnetic_flux_density)
# Note that this also means the vector of the (total) magnetic flux density is oriented in the
# positive direction of the `z`-axis.


@validate_input(relative_permeability_=absolute_permeability,
    current_=current,
    first_angle_=first_angle,
    second_angle_=second_angle,
    distance_=radial_distance)
@validate_output(magnetic_flux_density)
def calculate_induction(
    relative_permeability_: float,
    current_: Quantity,
    first_angle_: float | Quantity,
    second_angle_: float | Quantity,
    distance_: Quantity,
) -> Quantity:
    result_expr = solve(law, magnetic_flux_density, dict=True)[0][magnetic_flux_density]
    result_expr = result_expr.subs({
        absolute_permeability: relative_permeability_,
        current: current_,
        first_angle: first_angle_,
        second_angle: second_angle_,
        radial_distance: distance_
    })
    return Quantity(result_expr)

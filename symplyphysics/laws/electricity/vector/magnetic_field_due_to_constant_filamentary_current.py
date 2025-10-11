"""
Magnetic field due to constant filamentary current
==================================================

Known as the **Biot—Savart law**, it is an equation describing the magnetic field due to a constant
electric current.

**Notation:**

#. :quantity_notation:`vacuum_permeability`.

**Notes:**

#. This version of the law deals with a current in an infinitely thin wire. For a conductor of a
   finite thickness, the following relation must be used:

   .. math::

       I d \\vec{\\ell} = \\vec{J} dV

#. To find the total magnetic flux density, calculate the line integral over the whole contour.

**Conditions:**

#. The system is in a vacuum.

**Links:**

#. `Wikipedia — Biot—Savart law <https://en.wikipedia.org/wiki/Biot%E2%80%93Savart_law>`__.
"""

from sympy import Eq, pi
from symplyphysics import Quantity, validate_input, validate_output, symbols

from symplyphysics.core.experimental.vectors import clone_as_vector_symbol, VectorCross, VectorNorm
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

magnetic_flux_density_change = clone_as_vector_symbol(
    symbols.magnetic_flux_density,
    display_symbol="dB",
    display_latex="d \\vec{B}",
)
"""
Infinitesimal change of the :symbols:`magnetic_flux_density` at a given point in space.
"""

absolute_permeability = symbols.absolute_permeability
"""
:symbols:`absolute_permeability` of the medium.
"""

current = symbols.current
"""
:symbols:`current` in the contour.
"""

position_vector = clone_as_vector_symbol(symbols.distance_to_origin)
"""
Position vector of the point at which the magnetic flux density is measured. Also see
:symbols:`distance_to_origin`.
"""

contour_element_position_vector = clone_as_vector_symbol(symbols.distance_to_origin,
    display_symbol="l",
    display_latex="\\vec{\\ell}")
"""
Position vector of a point on the integration path. Also see :symbols:`distance_to_origin`.
"""

# FIXME: Better use `symplyphysics.core.operations.symbolic.ExactDifferential`, but this raises an
# error that it isn't a vector within VectorCross and VectorNorm.
contour_element_displacement = clone_as_vector_symbol(
    symbols.euclidean_distance,
    display_symbol="dl",
    display_latex="d \\vec{\\ell}",
)
"""
A vector along the integration path whose magnitude is the length of the differential element in
the direction of conventional current. Also see :symbols:`euclidean_distance`.
"""

law = Eq(
    magnetic_flux_density_change,
    (absolute_permeability / (4 * pi)) * ((current *
    VectorCross(contour_element_displacement, position_vector - contour_element_position_vector)) /
    VectorNorm(position_vector - contour_element_position_vector)**3),
)
"""
:laws:symbol::

:laws:latex::
"""

# Derive from Maxwell's equations and stationarity conditions
# Link in Russian: <https://ru.wikipedia.org/wiki/%D0%97%D0%B0%D0%BA%D0%BE%D0%BD_%D0%91%D0%B8%D0%BE_%E2%80%94_%D0%A1%D0%B0%D0%B2%D0%B0%D1%80%D0%B0_%E2%80%94_%D0%9B%D0%B0%D0%BF%D0%BB%D0%B0%D1%81%D0%B0#%D0%92%D1%8B%D0%B2%D0%BE%D0%B4_%D0%B7%D0%B0%D0%BA%D0%BE%D0%BD%D0%B0_%D0%B8%D0%B7_%D1%83%D1%80%D0%B0%D0%B2%D0%BD%D0%B5%D0%BD%D0%B8%D0%B9_%D0%9C%D0%B0%D0%BA%D1%81%D0%B2%D0%B5%D0%BB%D0%BB%D0%B0>


@validate_input(
    absolute_permeability_=absolute_permeability,
    current_=current,
    position_vector_=position_vector,
    contour_element_position_vector_=contour_element_position_vector,
    contour_element_displacement_=contour_element_displacement,
)
@validate_output(magnetic_flux_density_change)
def calculate_magnetic_flux_density_change(
    absolute_permeability_: Quantity,
    current_: Quantity,
    position_vector_: QuantityCoordinateVector,
    contour_element_position_vector_: QuantityCoordinateVector,
    contour_element_displacement_: QuantityCoordinateVector,
) -> QuantityCoordinateVector:
    result = law.rhs.subs({
        absolute_permeability: absolute_permeability_,
        current: current_,
        position_vector: position_vector_,
        contour_element_position_vector: contour_element_position_vector_,
        contour_element_displacement: contour_element_displacement_,
    })

    for norm in result.atoms(VectorNorm):
        arg = QuantityCoordinateVector.from_expr(norm.args[0])
        result = result.subs(norm, VectorNorm(arg))

    return QuantityCoordinateVector.from_expr(result)

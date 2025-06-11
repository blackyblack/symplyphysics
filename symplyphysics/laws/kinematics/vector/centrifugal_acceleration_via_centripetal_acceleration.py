"""
Centrifugal acceleration via centripetal acceleration
=====================================================

The vector of *centrifugal acceleration* has the same magnitude as the vector of *centripetal
acceleration* but is directed oppositely to it.

**Links:**

#. `BYJU's <https://byjus.com/physics/centripetal-and-centrifugal-force/>`__.
"""

from sympy import Eq
from symplyphysics import validate_input, validate_output, symbols

from symplyphysics.core.experimental.vectors import clone_as_vector_symbol
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

centrifugal_acceleration = clone_as_vector_symbol(
    symbols.acceleration,
    display_symbol="a_cf",
    display_latex="{\\vec a}_\\text{cf}",
)
"""
Vector of centrifugal :symbols:`acceleration` experienced by the body in the non-inertial,
rotating frame.
"""

centripetal_acceleration = clone_as_vector_symbol(
    symbols.acceleration,
    display_symbol="a_cp",
    display_latex="{\\vec a}_\\text{cp}",
)
"""
Vector of centripetal :symbols:`acceleration` experienced by the rotating body in the inertial
frame.
"""

law = Eq(centrifugal_acceleration, -1 * centripetal_acceleration)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(centripetal_acceleration_=centripetal_acceleration)
@validate_output(centrifugal_acceleration)
def calculate_centrifugal_acceleration(
        centripetal_acceleration_: QuantityCoordinateVector) -> QuantityCoordinateVector:
    result = law.rhs.subs(centripetal_acceleration, centripetal_acceleration_)

    return QuantityCoordinateVector.from_expr(result)

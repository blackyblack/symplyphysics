"""
Vector area is unit normal times scalar area
============================================

**Vector area**, or **oriented area**, is a vector quantity equal to the surface integral of the
surface normal. If the unit normal is constant at all points of the surface, integral can be
reduced to a product of the unit normal to the scalar area of the surface.

**Notes:**

#. If the normal changes direction across the surface, you can divide the surface into parts of
   constant unit normal and sum up the vector areas of all parts:

   .. math::

       \\vec A = \\sum_i \\vec A_i = \\sum_i \\vec n_i A_i

   Alternatively, use the surface integral that sums up infinitesimal vector areas across the
   surface:

   .. math::

       \\vec A = \\iint \\limits_S d \\vec A = \\iint \\limits_S \\vec n (\\vec r) dA

**Conditions:**

#. The surface is bounded (i.e. finite).

#. The surface normal is the same thoughout the given region, which can be achieved by choosing a
   small enough (or infinitesimal) surface.

**Links:**

#. `Wikipedia — Vector area <https://en.wikipedia.org/wiki/Vector_area#Definition>`__.
"""

from typing import Optional

from sympy import Eq
from symplyphysics import symbols, Quantity, validate_input, validate_output, assert_equal

from symplyphysics.core.experimental.vectors import clone_as_vector_symbol, VectorSymbol, VectorNorm
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

vector_area = clone_as_vector_symbol(symbols.area)
"""
Vector :symbols:`area` pertaining to the given region.
"""

unit_normal = VectorSymbol("n")
"""
Unit vector `normal <https://en.wikipedia.org/wiki/Normal_(geometry)#>`__ to the surface.

**Notes:**

#. :math:`\\left \\Vert \\vec n \\right \\Vert = 1`, i.e. it has unit magnitude.
"""

scalar_area = symbols.area
"""
Scalar :symbols:`area` of the given region.
"""

law = Eq(vector_area, unit_normal * scalar_area)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    unit_normal_=unit_normal,
    scalar_area_=scalar_area,
)
@validate_output(vector_area)
def calculate_vector_area(
    unit_normal_: QuantityCoordinateVector,
    scalar_area_: Quantity,
    relative_tolerance_: Optional[float] = None,
) -> QuantityCoordinateVector:
    assert_equal(VectorNorm(unit_normal_), 1, relative_tolerance=relative_tolerance_)

    result = law.rhs.subs({
        unit_normal: unit_normal_,
        scalar_area: scalar_area_,
    })

    return QuantityCoordinateVector.from_expr(result)

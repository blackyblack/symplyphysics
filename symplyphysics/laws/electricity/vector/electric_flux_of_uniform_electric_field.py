"""
Electric flux of uniform electric field
=======================================

Electric field at a point in space can be found by placing there a test charge and measuring
the electrostatic force that is applied to it.

**Notation:**

#. :math:`\\left( \\vec a, \\vec b \\right)` (:code:`dot(a, b)`) is the dot product between vectors
   :math:`\\vec a` and :math:`\\vec b`.

**Notes:**

#. Vector area is a vector quantity whose magnitude denotes the area of the surface it
   represents and whose direction denotes the orientation of the surface.

**Conditions:**

#. The electric field is uniform. This might be achieved by choosing a small enough surface
   that the electric field would be constant throughout it.

**Links:**

#. `Electric flux <https://en.wikipedia.org/wiki/Electric_flux#Overview>`__.
"""

from sympy import Eq
from symplyphysics import Quantity, validate_input, validate_output, symbols

from symplyphysics.core.experimental.vectors import clone_as_vector_symbol, VectorDot
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

electric_flux = symbols.electric_flux
"""
:symbols:`electric_flux`.
"""

electric_field = clone_as_vector_symbol(symbols.electric_field_strength)
"""
Vector of the electric field. See :symbols:`electric_field_strength`.
"""

area = clone_as_vector_symbol(symbols.area)
"""
Area pseudovector, i.e. a vector that is aligned in the direction of the unit normal to the surface
and whose magnitude is equal to the area of the surface. See :symbols:`area`.
"""

law = Eq(electric_flux, VectorDot(electric_field, area))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(electric_field_=electric_field, area_=area)
@validate_output(electric_flux)
def calculate_electric_flux(
    electric_field_: QuantityCoordinateVector,
    area_: QuantityCoordinateVector,
) -> Quantity:
    result = law.rhs.subs({
        electric_field: electric_field_,
        area: area_,
    })

    return Quantity(result)


# UNIQUE_LAW_ID: 527

r"""
Electric flux of uniform electric field
=======================================

Electric field at a point in space can be found by placing there a test charge and measuring
the electrostatic force that is applied to it.

**Notation:**

#. :math:`\vec a \cdot \vec b` (:code:`dot(a, b)`) is the dot product between vectors
   :math:`\vec a` and :math:`\vec b`.

**Notes:**

#. Vector area is a vector quantity whose magnitude denotes the area of the surface it
   represents and whose direction denotes the orientation of the surface.

**Conditions:**

#. The electric field is uniform. This might be achieved by choosing a small enough surface
   that the electric field would be constant throughout it.

**Links:**

#. `Electric flux <https://en.wikipedia.org/wiki/Electric_flux#Overview>`__.
"""

from sympy import Expr
from symplyphysics import (
    Quantity,
    QuantityVector,
    dot_vectors,
    units,
    validate_input,
    validate_output,
    Vector,
)

def electric_flux_law(
    electric_field_: Vector,
    area_: Vector,
) -> Expr:
    r"""
    Electric flux via electric field and vector area.

    Law:
        :code:`Phi_E = dot(E, A)`

    Latex:
        .. math::
            \Phi_E = \vec E \cdot \vec A

    :param electric\_field\_: vector of electric field

        Symbol: :code:`E`

        Latex: :math:`\vec E`

        Dimension: :code:`voltage/length`

    :param area\_: vector area

        Symbol: :code:`A`

        Latex: :math:`\vec A`

        Dimension: :code:`area`

    :return: electric flux

        Symbol: :code:`Phi_E`

        Latex: :math:`\Phi_E`

        Dimension: :code:`voltage*length`
    """

    return dot_vectors(electric_field_, area_)


@validate_input(
    electric_field_=units.voltage / units.length,
    area_=units.area,
)
@validate_output(units.voltage * units.length)
def calculate_electric_flux(
    electric_field_: QuantityVector,
    area_: QuantityVector,
) -> Quantity:
    result = electric_flux_law(
        electric_field_=electric_field_.to_base_vector(),
        area_=area_.to_base_vector(),
    )
    return Quantity(result)

"""
Hydraulic stress is bulk modulus times strain
=============================================

When an object undergoes hydraulic compression due to a stress exerted by a surrounding liquid,
the pressure (hydraulic stress) on the object due to the fluid is proportional to the fractional
change in the object's volume due to that pressure and the bulk modulus of the object. Thus, bulk
modulus of a substance is a measure of its resistance to bulk compression.
"""

from sympy import Eq
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

hydraulic_stress = Symbol("hydraulic_stress", units.pressure)
"""
Hydraulic stress.

Symbol:
    :code:`p`
"""

bulk_modulus = Symbol("bulk_modulus", units.pressure)
"""
Bulk modulus of the material.

Symbol:
    :code:`B`
"""

fractional_volume_change = Symbol("fractional_volume_change", dimensionless)
r"""
:doc:`Fractional volume change <laws.quantities.fractional_change_is_change_over_initial_value>`.

Symbol:
    :code:`e_V`

Latex:
    :math:`e_V`
"""

law = Eq(hydraulic_stress, bulk_modulus * fractional_volume_change)
r"""
:code:`p = B * e_V`

Latex:
    .. math::
        p = B e_V
"""


@validate_input(
    bulk_modulus_=bulk_modulus,
    fractional_volume_change_=fractional_volume_change,
)
@validate_output(hydraulic_stress)
def calculate_hydraulic_stress(
    bulk_modulus_: Quantity,
    fractional_volume_change_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        bulk_modulus: bulk_modulus_,
        fractional_volume_change: fractional_volume_change_,
    })
    return Quantity(result)

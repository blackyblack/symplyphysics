"""
Hydraulic stress is bulk modulus times strain
=============================================

When an object undergoes hydraulic compression due to a stress exerted by a surrounding liquid,
the pressure (hydraulic stress) on the object due to the fluid is proportional to the fractional
change in the object's volume due to that pressure and the bulk modulus of the object. Thus, bulk
modulus of a substance is a measure of its resistance to bulk compression.

**Links:**

#. `Wikipedia, derivable from 12.5.3 <https://phys.libretexts.org/Workbench/PH_245_Textbook_V2/12%3A_Static_Equilibrium_and_Elasticity/12.05%3A_Stress_Strain_and_Elastic_Modulus_(Part_2)>`__.

..
    TODO: rename `hydraulic stress` to `bulk stress`
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

bulk_stress = clone_as_symbol(symbols.pressure, display_symbol="Delta(p)", display_latex="\\Delta p")
"""
Bulk stress. See :symbols:`pressure`.
"""

bulk_modulus = symbols.bulk_modulus
"""
:symbols:`bulk_modulus` of the material.
"""

fractional_volume_change = clone_as_symbol(symbols.fractional_change, subscript="V")
"""
:symbols:`fractional_change` of volume. See :ref:`Fractional change is change over
initial value`.
"""

law = Eq(bulk_stress, bulk_modulus * fractional_volume_change)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    bulk_modulus_=bulk_modulus,
    fractional_volume_change_=fractional_volume_change,
)
@validate_output(bulk_stress)
def calculate_hydraulic_stress(
    bulk_modulus_: Quantity,
    fractional_volume_change_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        bulk_modulus: bulk_modulus_,
        fractional_volume_change: fractional_volume_change_,
    })
    return Quantity(result)

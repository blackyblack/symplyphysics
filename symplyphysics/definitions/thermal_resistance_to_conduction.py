"""
Thermal resistance to conduction
================================

Thermal resistance to conduction is a physical quantity that was introduced in the engineering
practice for insulators: greater values of thermal resistance mean better insulation properties of
a material of given thickness.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

thermal_insulance = symbols.thermal_insulance
"""
:symbols:`thermal_insulance` of the insulator, also called an **R-value**.
"""

slab_thickness = symbols.thickness
"""
:symbols:`thickness` of the insulator.
"""

thermal_conductivity = symbols.thermal_conductivity
"""
:symbols:`thermal_conductivity` of the insulating material.
"""

definition = Eq(thermal_insulance, slab_thickness / thermal_conductivity)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    slab_thickness_=slab_thickness,
    thermal_conductivity_=thermal_conductivity,
)
@validate_output(thermal_insulance)
def calculate_thermal_insulance(
    slab_thickness_: Quantity,
    thermal_conductivity_: Quantity,
) -> Quantity:
    result = definition.rhs.subs({
        slab_thickness: slab_thickness_,
        thermal_conductivity: thermal_conductivity_,
    })
    return Quantity(result)

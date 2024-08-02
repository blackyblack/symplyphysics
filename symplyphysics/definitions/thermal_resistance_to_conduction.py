"""
Thermal resistance to conduction
================================

Thermal resistance to conduction is a physical quantity that was introduced in the engineering
practice for insulators: greater values of thermal resistance mean better insulation properties of
a material of given thickness.
"""

from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

# Definition: R = L/k
## R - thermal resistance to conduction, also known as R-value
## L - thickness of slab
## k - thermal conductivity of material

thermal_resistance = Symbol("thermal_resistance", units.area * units.temperature / units.power)
"""
Thermal resistance to conduction of the insulator, also called R-value.

Symbol:
    :code:`R`
"""

slab_thickness = Symbol("slab_thickness", units.length)
"""
Thickness of the insulator.

Symbol:
    :code:`l`
"""

thermal_conductivity = Symbol("thermal_conductivity",
    units.power / (units.length * units.temperature))
"""
Thermal conductivity of the insulating material.

Symbol:
    :code:`k`
"""

definition = Eq(thermal_resistance, slab_thickness / thermal_conductivity)
r"""
:code:`R = l / k`

Latex:
    .. math::
        R = \frac{l}{k}
"""


@validate_input(
    slab_thickness_=slab_thickness,
    thermal_conductivity_=thermal_conductivity,
)
@validate_output(thermal_resistance)
def calculate_thermal_resistance(
    slab_thickness_: Quantity,
    thermal_conductivity_: Quantity,
) -> Quantity:
    result = definition.rhs.subs({
        slab_thickness: slab_thickness_,
        thermal_conductivity: thermal_conductivity_,
    })
    return Quantity(result)

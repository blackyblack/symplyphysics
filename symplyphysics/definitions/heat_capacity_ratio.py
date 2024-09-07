"""
Heat capacity ratio
===================

The *heat capacity ratio*, also known as the *adiabatic index*, the *ratio of specific heats*, or
the *insentropic expansion factor*, is the ratio of the heat capacity at constant pressure
to that of constant volume. The heat capacity ratio is used in the description of thermodynamic
reversible processes; the speed of sound also depends on this factor.

**Notes:**

#. One can also use intensive heat capacities, such as specific or molar ones, in place of the 
   extensive heat capacity presented here.
"""

from sympy import Eq
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
    convert_to_float,
)

heat_capacity_ratio = SymbolNew("gamma", dimensionless, display_latex="\\gamma")
"""
Heat capacity ratio of the system.
"""

isobaric_heat_capacity = SymbolNew("C_p", units.energy / units.temperature)
"""
Heat capacity of the system at constant pressure.
"""

isochoric_heat_capacity = SymbolNew("C_V", units.energy / units.temperature)
"""
Heat capacity of the system at constant volume.
"""

definition = Eq(heat_capacity_ratio, isobaric_heat_capacity / isochoric_heat_capacity)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    isobaric_heat_capacity_=isobaric_heat_capacity,
    isochoric_heat_capacity_=isochoric_heat_capacity,
)
@validate_output(heat_capacity_ratio)
def calculate_heat_capacity_ratio(
    isobaric_heat_capacity_: Quantity,
    isochoric_heat_capacity_: Quantity,
) -> float:
    result = definition.rhs.subs({
        isobaric_heat_capacity: isobaric_heat_capacity_,
        isochoric_heat_capacity: isochoric_heat_capacity_,
    })
    return convert_to_float(result)

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
    Symbol,
    validate_input,
    validate_output,
    convert_to_float,
)


heat_capacity_ratio = Symbol("heat_capacity_ratio", dimensionless)
r"""
Heat capacity ratio of the system.

Symbol:
    :code:`gamma`

Latex:
    :math:`\gamma`
"""

isobaric_heat_capacity = Symbol("isobaric_heat_capacity", units.energy / units.temperature)
"""
Heat capacity of the system at constant pressure.

Symbol:
    :code:`C_p`

Latex:
    :math:`C_p`
"""

isochoric_heat_capacity = Symbol("isochoric_heat_capacity", units.energy / units.temperature)
"""
Heat capacity of the system at constant volume.

Symbol:
    :code:`C_V`

Latex:
    :math:`C_V`
"""

definition = Eq(heat_capacity_ratio, isobaric_heat_capacity / isochoric_heat_capacity)
r"""
:code:`gamma = C_p / C_V`

Latex:
    .. math::
        \gamma = \frac{C_p}{C_V}
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

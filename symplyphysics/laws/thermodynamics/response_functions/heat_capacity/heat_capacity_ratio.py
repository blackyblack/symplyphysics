"""
Heat capacity ratio
===================

The *heat capacity ratio* (also called the *adiabatic index*, *ratio of specific heats*, or *isentropic expansion factor*)
is the ratio of heat capacity at constant pressure to heat capacity at constant volume. It governs adiabatic processes and
influences the speed of sound in a medium.

**Notes:**

#. One can also use intensive heat capacities, e.g. specific or molar, in place of the extensive heat capacities presented here.

**Links:**

#. `Wikipedia – Heat capacity ratio <https://en.wikipedia.org/wiki/Heat_capacity_ratio>`__
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    convert_to_float,
    symbols,
    clone_as_symbol,
)

heat_capacity_ratio = symbols.adiabatic_index
"""
:symbols:`adiabatic_index` of the system.
"""

isobaric_heat_capacity = clone_as_symbol(symbols.heat_capacity,
    display_symbol="C_p",
    display_latex="C_p")
"""
:symbols:`heat_capacity` of the system at constant pressure.
"""

isochoric_heat_capacity = clone_as_symbol(symbols.heat_capacity,
    display_symbol="C_V",
    display_latex="C_V")
"""
:symbols:`heat_capacity` of the system at constant volume.
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

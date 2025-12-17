"""
Isochoric and isobaric heat capacities of homogeneous substance
===============================================================

The **Mayer's relation** is the relation between heat capacity at constant pressure and heat
capacity at constant volume. In the current form it is applicable to any homogeneous substance,
not just ideal gases.

**Links:**

#. `Wikipedia, second formula <https://en.wikipedia.org/wiki/Mayer%27s_relation>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

isobaric_heat_capacity = clone_as_symbol(
    symbols.heat_capacity,
    display_symbol="C_p",
    display_latex="C_p",
)
"""
:symbols:`heat_capacity` at constant :symbols:`pressure`.
"""

isochoric_heat_capacity = clone_as_symbol(
    symbols.heat_capacity,
    display_symbol="C_V",
    display_latex="C_V",
)
"""
:symbols:`heat_capacity` at constant :symbols:`volume`.
"""

volume = symbols.volume
"""
:symbols:`volume` of the substance.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the substance.
"""

thermal_expansion_coefficient = clone_as_symbol(
    symbols.thermal_expansion_coefficient,
    subscript="V",
)
"""
:symbols:`thermal_expansion_coefficient` of the substance. Also see :doc:`Thermal volumetric expansion
coefficient <definitions.volumetric_coefficient_of_thermal_expansion>`.
"""

isothermal_compressibility = clone_as_symbol(
    symbols.thermodynamic_compressibility,
    subscript="T",
)
"""
:symbols:`thermodynamic_compressibility` of the substance. Also see :doc:`Isothermal compressibility
<definitions.thermodynamic_compressibility>`.
"""

law = Eq(
    isobaric_heat_capacity - isochoric_heat_capacity,
    volume * temperature * thermal_expansion_coefficient**2 / isothermal_compressibility,
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    isobaric_heat_capacity_=isobaric_heat_capacity,
    volume_=volume,
    temperature_=temperature,
    thermal_expansion_coefficient_=thermal_expansion_coefficient,
    isothermal_compressibility_=isothermal_compressibility,
)
@validate_output(isochoric_heat_capacity)
def calculate_isochoric_heat_capacity(
    isobaric_heat_capacity_: Quantity,
    volume_: Quantity,
    temperature_: Quantity,
    thermal_expansion_coefficient_: Quantity,
    isothermal_compressibility_: Quantity,
) -> Quantity:
    expr = solve(law, isochoric_heat_capacity)[0]
    result = expr.subs({
        isobaric_heat_capacity: isobaric_heat_capacity_,
        volume: volume_,
        temperature: temperature_,
        thermal_expansion_coefficient: thermal_expansion_coefficient_,
        isothermal_compressibility: isothermal_compressibility_,
    })
    return Quantity(result)

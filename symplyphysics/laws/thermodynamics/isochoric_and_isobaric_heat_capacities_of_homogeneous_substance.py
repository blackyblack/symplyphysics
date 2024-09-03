"""
Isochoric and isobaric heat capacities of homogeneous substance
===============================================================

The *Mayer's relation* is the relation between heat capacity at constant pressure and heat
capacity at constant volume. In the current form it is applicable to any homogeneous substance,
not just ideal gases.
"""

from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
)

isobaric_heat_capacity = Symbol("isobaric_heat_capacity", units.energy / units.temperature)
r"""
Heat capacity at constant pressure.

Symbol:
    :code:`C_p`

Latex:
    :math:`C_p`
"""

isochoric_heat_capacity = Symbol("isochoric_heat_capacity", units.energy / units.temperature)
r"""
Heat capacity at constant volume.

Symbol:
    :code:`C_V`

Latex:
    :math:`C_V`
"""

volume = Symbol("volume", units.volume)
"""
Volume of the substance.

Symbol:
    :code:`V`
"""

temperature = symbols.thermodynamics.temperature
"""
Temperature of the substance.
"""

thermal_expansion_coefficient = Symbol("thermal_expansion_coefficient", 1 / units.temperature)
r"""
:doc:`Thermal volumetric expansion coefficient <definitions.volumetric_coefficient_of_thermal_expansion>`
of the substance.

Symbol:
    :code:`alpha_V`

Latex:
    :math:`\alpha_V`
"""

isothermal_compressibility = Symbol("isothermal_compressibility", 1 / units.pressure)
r"""
:doc:`Isothermal compressibility <definitions.thermodynamic_compressibility>` of the substance.

Symbol:
    :code:`beta_T`

Latex:
    :math:`\beta_T`
"""

law = Eq(
    isobaric_heat_capacity - isochoric_heat_capacity,
    volume * temperature * thermal_expansion_coefficient**2 / isothermal_compressibility,
)
r"""
:code:`C_p - C_V = V * T * (alpha_V)^2 / beta_T`

Latex:
    .. math::
        C_p - C_V = V T \frac{\alpha_V^2}{\beta_T}
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

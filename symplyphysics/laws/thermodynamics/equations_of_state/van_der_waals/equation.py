"""
Van der Waals equation
======================

To more accurately describe the behavior of real gases at low temperatures,
a Van der Waals gas model was created, taking into account the forces of intermolecular interaction.
In this model, internal energy becomes a function not only of temperature, but also of molar volume.

The Van der Waals equation is one of the well-known approximate equations of state describing
the properties of a real gas, having a compact form and taking
into account the main characteristics of a gas with intermolecular interaction.

**Notation:**

#. :math:`R` is the molar gas constant.
"""

from sympy import (Eq, solve)
from symplyphysics import (symbols, units, Quantity, Symbol, validate_input,
    validate_output)

pressure = Symbol("pressure", units.pressure)
"""
Pressure in the van der Waals fluid.

Symbol:
    :code:`p`
"""

molar_volume = Symbol("molar_volume", units.volume / units.amount_of_substance)
"""
Volume of the van der Waals fluid per unit amount of substance.

Symbol:
    :code:`V_m`

Latex:
    :math:`V_m`
"""

temperature = symbols.thermodynamics.temperature
"""
:attr:`~symplyphysics.symbols.thermodynamics.temperature` of the van der Waals fluid.

Symbol:
    :code:`T`
"""

attractive_forces_parameter = Symbol(
    "attractive_forces_parameter",
    units.pressure * (units.volume / units.amount_of_substance)**2,
)
"""
Parameter of the van der Waals equation denoting the magnitude of attractive
forces between gas molecules.

Symbol:
    :code:`a`
"""

excluded_volume_parameter = Symbol(
    "excluded_volume_parameter",
    units.volume / units.amount_of_substance,
)
"""
Parameter of the van der Waals equation denoting an excluded molar volume
due to a finite size of molecules.

Symbol:
    :code:`b`
"""

law = Eq(
    (pressure + attractive_forces_parameter / molar_volume**2) * (molar_volume - excluded_volume_parameter),
    units.molar_gas_constant * temperature,
)
r"""
:code:`(p + a / V_m^2) * (V_m - b) = R * T`

Latex:
    .. math::
        \left( p + \frac{a}{V_m^2} \right) (V_m - b) = R T
"""


@validate_input(
    molar_volume_=molar_volume,
    temperature_=temperature,
    bonding_forces_parameter_=attractive_forces_parameter,
    molecules_volume_parameter_=excluded_volume_parameter,
)
@validate_output(pressure)
def calculate_pressure(
    molar_volume_: Quantity,
    temperature_: Quantity,
    bonding_forces_parameter_: Quantity,
    molecules_volume_parameter_: Quantity,
) -> Quantity:
    solved = solve(law, pressure, dict=True)[0][pressure]
    result_expr = solved.subs({
        molar_volume: molar_volume_,
        temperature: temperature_,
        attractive_forces_parameter: bonding_forces_parameter_,
        excluded_volume_parameter: molecules_volume_parameter_
    })
    return Quantity(result_expr)

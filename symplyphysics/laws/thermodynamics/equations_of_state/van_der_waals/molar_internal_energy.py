"""
Molar internal energy
=====================

If the equation of state is known, the internal energy of a substance can be found
as a function of volume at constant temperature.

**Conditions:**

#. The fluid is homogeneous and in a single phase state.
"""

from sympy import Eq, Integral
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    Function,
    validate_input,
    validate_output,
    symbols,
)

molar_internal_energy = Symbol("molar_internal_energy", units.energy / units.amount_of_substance)
"""
Internal energy of the van der Waals fluid per unit amount of substance.

Symbol:
    :code:`u`
"""

isochoric_molar_heat_capacity = Function(
    "isochoric_molar_heat_capacity",
    units.energy / (units.temperature * units.amount_of_substance),
)
r"""
Heat capacity at constant volume per unit amount of substance.

Symbol:
    :code:`c_V(T)`

Latex:
    :math:`c_V(T)`
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

molar_volume = Symbol("molar_volume", units.volume / units.amount_of_substance)
r"""
Volume of the van der Waals fluid per unit amount of substance.

Symbol:
    :code:`V_m`

Latex:
    :math:`V_m`
"""

law = Eq(
    molar_internal_energy,
    Integral(isochoric_molar_heat_capacity(temperature), temperature) -
    attractive_forces_parameter / molar_volume,
)
r"""
:code:`u = Integral(c_V(T), T) - a / V_m`

Latex:
    .. math::
        u = \int c_V(T) \, dT - \frac{a}{V_m}
"""


@validate_input(
    isochoric_molar_heat_capacity_=isochoric_molar_heat_capacity,
    temperature_=temperature,
    bonding_forces_parameter_=attractive_forces_parameter,
    molar_volume_=molar_volume,
)
@validate_output(molar_internal_energy)
def calculate_internal_energy(
    isochoric_molar_heat_capacity_: Quantity,
    temperature_: Quantity,
    bonding_forces_parameter_: Quantity,
    molar_volume_: Quantity,
) -> Quantity:
    # Note that internal energy is only known up to a constant term
    # Isochoric heat capacity is assumed to be a constant independent of temperature

    isochoric_molar_heat_capacity_function = isochoric_molar_heat_capacity_
    result = law.rhs.subs(isochoric_molar_heat_capacity(temperature),
        isochoric_molar_heat_capacity_function).doit().subs({
        temperature: temperature_,
        attractive_forces_parameter: bonding_forces_parameter_,
        molar_volume: molar_volume_,
        })
    return Quantity(result)

r"""
Dieterici equation
==================

*Dieterici equation* is another type of semi-empirical equations approximating real gases
along with the more well-known van der Waals equation of state.

**Notation:**

#. :math:`R` is the molar gas constant.

**Notes:**

#. Like the van der Waals equation of state, the Dieterici equation is semi-empirical.
#. It approximates moderate pressures of real gases much better than the van der Waals equation
   within the conditions stated below.
#. Can be converted to the van der Waals equation under an additional limit :math:`a \ll R T V_m`

**Conditions:**

#. Only applicable in the limits :math:`b \ll V_m` and :math:`a \ll p V_m^2`.
#. Inapplicable for large pressures.
"""

from sympy import Eq, solve, exp
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
)

pressure = Symbol("pressure", units.pressure)
"""
Pressure inside the system.

Symbol:
    :code:`p`
"""

molar_volume = Symbol("molar_volume", units.volume / units.amount_of_substance)
r"""
Volume of the system per amount of substance.

Symbol:
    :code:`V_m`

Latex:
    :math:`V_m`
"""

temperature = symbols.thermodynamics.temperature
"""
:attr:`~symplyphysics.symbols.thermodynamics.temperature` of the system.

Symbol:
    :code:`T`
"""

attractive_forces_parameter = Symbol(
    "attractive_forces_parameter",
    units.pressure * (units.volume / units.amount_of_substance)**2,
)
"""
Parameter specific to each individual substance, usually attributed to the magnitude of
attractive forces between particles of the system.

Symbol:
    :code:`a`
"""

excluded_volume_parameter = Symbol(
    "excluded_volume_parameter",
    units.volume / units.amount_of_substance,
)
"""
Parameter specific to each individual substance, usually attributed to the amount of
excluded molar volume due to a finite size of particles.
"""

law = Eq(
    pressure * (molar_volume - excluded_volume_parameter),
    units.molar_gas_constant * temperature * exp(-1 * attractive_forces_parameter /
    (units.molar_gas_constant * temperature * molar_volume)))
r"""
:code:`p * (V_m - b) = R * T * exp(-1 * a / (R * T * V_m))`

Latex:
    .. math::
        p \left( V_m - b \right) = R T \exp \left( - \frac{a}{R T V_m} \right)
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
    expr = solve(law, pressure)[0]
    result = expr.subs({
        molar_volume: molar_volume_,
        temperature: temperature_,
        attractive_forces_parameter: bonding_forces_parameter_,
        excluded_volume_parameter: molecules_volume_parameter_,
    })
    return Quantity(result)

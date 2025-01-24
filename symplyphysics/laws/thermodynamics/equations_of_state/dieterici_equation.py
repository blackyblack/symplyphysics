r"""
Dieterici equation
==================

*Dieterici equation* is another type of semi-empirical equations approximating real gases
along with the more well-known van der Waals equation of state.

**Notation:**

#. :quantity_notation:`molar_gas_constant`.

**Notes:**

#. Like the van der Waals equation of state, the Dieterici equation is semi-empirical.
#. It approximates moderate pressures of real gases much better than the van der Waals equation
   within the conditions stated below.
#. Can be converted to the van der Waals equation under an additional limit :math:`a \ll R T V_m`.

**Conditions:**

#. Only applicable in the limits :math:`b \ll V_m` and :math:`a \ll p V_m^2`. Refer to symbols below.
#. Inapplicable for high pressures.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Real_gas#Dieterici_model>`__.
"""

from sympy import Eq, solve, exp
from symplyphysics import (
    units,
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
    symbols,
    quantities,
)

pressure = symbols.pressure
"""
:symbols:`pressure` inside the system.
"""

molar_volume = SymbolNew("V_m", units.volume / units.amount_of_substance)
"""
:symbols:`volume` of the system per :symbols:`amount_of_substance`.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the system.
"""

attractive_forces_parameter = symbols.attractive_forces_parameter
"""
:symbols:`attractive_forces_parameter`.
"""

excluded_volume_parameter = symbols.excluded_volume_parameter
"""
:symbols:`excluded_volume_parameter`.
"""

law = Eq(
    pressure * (molar_volume - excluded_volume_parameter),
    quantities.molar_gas_constant * temperature * exp(-1 * attractive_forces_parameter /
    (quantities.molar_gas_constant * temperature * molar_volume)))
"""
:laws:symbol::

:laws:latex::
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

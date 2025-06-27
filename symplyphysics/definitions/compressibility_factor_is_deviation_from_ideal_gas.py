"""
Compressibility factor is deviation from ideal gas
==================================================

The **compressibility factor** (also called the **compression factor** or **gas deviation factor**) quantifies how far a real gas departs from ideal-gas behaviour.

**Notation:**

#. :quantity_notation:`molar_gas_constant`.

**Conditions:**

#. Temperature is strictly greater than zero.
#. Amount of substance, pressure, and volume are all finite and non-zero.
#. The gas sample is at thermodynamic equilibrium.

**Notes:**

#. Can be equivalently defined as the ratio of the molar volume :math:`V/n` of the real gas to the molar volume :math:`RT/\\rho` of the corresponding ideal gas at the same temperature and pressure.
#. :math:`Z = 1` corresponds to ideal-gas behaviour.
#. At high pressures repulsive interactions dominate, giving :math:`Z > 1`.
#. At low pressures attractive interactions dominate, giving :math:`Z < 1`.

**Links:**

#. `Wikipedia – Compressibility factor <https://en.wikipedia.org/wiki/Compressibility_factor#Definition_and_physical_significance>`__
"""

from sympy import Eq
from symplyphysics import (
    quantities,
    Quantity,
    validate_input,
    validate_output,
    symbols,
    convert_to_float,
)

compressibility_factor = symbols.compressibility_factor
"""
:symbols:`compressibility_factor` of the real gas.
"""

pressure = symbols.pressure
"""
:symbols:`pressure` of the gas.
"""

volume = symbols.volume
"""
:symbols:`volume` of the gas.
"""

amount_of_substance = symbols.amount_of_substance
"""
:symbols:`amount_of_substance` of the gas.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the gas.
"""

definition = Eq(compressibility_factor,
    (pressure * volume) / (amount_of_substance * quantities.molar_gas_constant * temperature))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    pressure_=pressure,
    volume_=volume,
    amount_of_substance_=amount_of_substance,
    temperature_=temperature,
)
@validate_output(compressibility_factor)
def calculate_compressibility_factor(
    pressure_: Quantity,
    volume_: Quantity,
    amount_of_substance_: Quantity,
    temperature_: Quantity,
) -> float:
    result = definition.rhs.subs({
        pressure: pressure_,
        volume: volume_,
        amount_of_substance: amount_of_substance_,
        temperature: temperature_,
    })
    return convert_to_float(result)

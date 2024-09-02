r"""
Compressibility factor is deviation from ideal gas
==================================================

The *compressibility factor*, also known as the *compression factor* or the *gas deviation factor*,
describes the deviation of a real gas from ideal gas behaviour. In general, the deviation from
ideal gas behaviour becomes more prominent the closer the gas is to a phase change, the lower
the temperature or the larger the pressure.

**Notation:**

#. :math:`R` is the molar gas constant.

**Notes:**

#. Can be equivalently defined as the ratio of the molar volume :math:`\frac{V}{n}` of the real gas to the
   molar volume :math:`\frac{R T}{\rho}` of the corresponding ideal gas at the same temperature and pressure.
#. :math:`Z = 1` is the case of ideal gas behaviour.
#. At high pressures molecules collide more often leading to an increase of repulsive forces between
   molecules, making the molar volume of the real gas greater than that of ideal gas, in other words the 
   particles have a larger extended volume, leading to :math:`Z > 1`.
#. At lower pressures, molecules are free to move and attractive forces dominate, leading to :math:`Z < 1`.
"""

from sympy import Eq
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
    convert_to_float,
)

compressibility_factor = Symbol("compressibility_factor", dimensionless)
"""
Compressibility factor of the real gas.

Symbol:
    :code:`Z`
"""

pressure = Symbol("pressure", units.pressure)
"""
Pressure of the gas.

Symbol:
    :code:`p`
"""

volume = Symbol("volume", units.volume)
"""
Volume of the gas.

Symbol:
    :code:`V`
"""

amount_of_substance = Symbol("amount_of_substance", units.amount_of_substance)
"""
Amount of gas substance.

Symbol:
    :code:`n`
"""

temperature = symbols.thermodynamics.temperature
"""
Gas :attr:`~symplyphysics.symbols.thermodynamics.temperature`.
"""

definition = Eq(compressibility_factor,
    (pressure * volume) / (amount_of_substance * units.molar_gas_constant * temperature))
r"""
:code:`Z = (p * V) / (n * R * T)`

Latex:
    .. math::
        Z = \frac{p V}{n R T}
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

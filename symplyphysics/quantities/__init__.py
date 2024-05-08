"""! @brief Contains useful physical quantities. Fundamental quantites can be also found in `sympy.physics.units` module."""

from sympy.physics import units
from ..core.symbols.quantities import Quantity

# Zero Celsius degrees. The temperature at which water freezes.
# It is also temperature for Standard Temperature and Pressure (STP)
standard_conditions_temperature = Quantity(273.15 * units.kelvin)

# Approximately 25 Celsius degrees. Commonly used temperature for tabulation
# purposes.
standard_laboratory_temperature = Quantity(298 * units.kelvin)

# Mass of stationary electron
electron_rest_mass = Quantity(9.1093837015e-31 * units.kilogram)

__all__ = [
    "standard_conditions_temperature",
]

"""! @brief Contains useful physical quantities. Fundamental quantites can be also found in `sympy.physics.units` module."""

from symplyphysics import Quantity, units


# Zero Celsius degrees. The temperature at which water freezes.
# It is also temperature for Standard Temperature and Pressure (STP)
zero_celsius = Quantity(273.15 * units.kelvin)

# Approximately 25 Celsius degrees. Commonly used temperature for tabulation
# purposes.
standard_conditions_temperature = Quantity(298 * units.kelvin)

__all__ = [
    "zero_celsius",
]

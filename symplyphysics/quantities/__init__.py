"""! @brief Contains useful physical quantities. Fundamental quantites can be also found in `sympy.physics.units` module."""

from sympy.physics import units
from ..core.symbols.quantities import Quantity

# Zero Celsius degrees. The temperature at which water freezes.
# It is also temperature for Standard Temperature and Pressure (STP)
zero_celsius = Quantity(273.15 * units.kelvin)

# Approximately 25 Celsius degrees. Commonly used temperature for tabulation
# purposes.
standard_conditions_temperature = Quantity(298 * units.kelvin)

# A physical constant commonly used in thermodynamics and statistical physics.
# It is the proportionality factor relating the average relative thermal energy
# of particles in a gas with the thermodynamic temperature of the gas.
boltzmann_constant = Quantity(1.380649e-23 * units.joule / units.kelvin)

__all__ = [
    "zero_celsius",
    "boltzmann_constant",
]

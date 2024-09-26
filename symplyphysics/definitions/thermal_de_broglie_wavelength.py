"""
Thermal de Broglie wavelength
=============================

The thermal de Broglie wavelength can be roughly described as the average de Broglie wavelength
of particles in an ideal gas at a specified temperature. When compared to average inter-particle
spacing in the gas, it can be used to tell if the gas can be considered to be a classical or
Maxwell-Boltzmann gas, in which case the thermal wavelength must be much smaller than the average
inter-particle spacing. Otherwise, quantum effects must be taken into account.

**Notation:**

#. :quantity_notation:`hbar`.
#. :quantity_notation:`boltzmann_constant`.
"""

from sympy import Eq, pi
from symplyphysics import (
    quantities,
    Quantity,
    validate_input,
    validate_output,
    symbols,
    sqrt,
)

thermal_wavelength = symbols.wavelength
"""
Thermal de Broglie :symbols:`wavelength` of the gas.
"""

mass = symbols.mass
"""
:symbols:`mass` of a single gas particle.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the gas.
"""

definition = Eq(thermal_wavelength,
    quantities.hbar * sqrt(2 * pi / (mass * quantities.boltzmann_constant * temperature)))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    mass_=mass,
    temperature_=temperature,
)
@validate_output(thermal_wavelength)
def calculate_thermal_wavelength(
    mass_: Quantity,
    temperature_: Quantity,
) -> Quantity:
    result = definition.rhs.subs({
        mass: mass_,
        temperature: temperature_,
    })
    return Quantity(result)

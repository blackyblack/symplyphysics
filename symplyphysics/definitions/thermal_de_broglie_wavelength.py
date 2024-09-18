r"""
Thermal de Broglie wavelength
=============================

The thermal de Broglie wavelength can be roughly described as the average de Broglie wavelength
of particles in an ideal gas at a specified temperature. When compared to average interparticle
spacing in the gas, it can be used to tell if the gas can be considered to be a classical or
Maxwell-Boltzmann gas, in which case the thermal wavelength must be much smaller than the average
interparticle spacing. Otherwise, quantum effects must be taken into account.

**Notation:**

#. :quantity_notation:`hbar`.
#. :quantity_notation:`boltzmann_constant`.
"""

from sympy import Eq, sqrt, pi
from symplyphysics import (
    quantities,
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
)

thermal_wavelength = Symbol("thermal_wavelength", units.length)
r"""
Thermal de Broglies wavelength of the gas.

Symbol:
    :code:`lambda`

Latex:
    :math:`\lambda`
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
r"""
:code:`lambda = hbar * sqrt(2 * pi / (m * k_B * T))`

Latex:
    .. math::
        \lambda = \hbar \sqrt{\frac{2 \pi}{m k_\text{B} T}}
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

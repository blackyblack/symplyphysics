r"""
Thermal de Broglie wavelength
=============================

The thermal de Broglie wavelength can be roughly described as the average de Broglie wavelength
of particles in an ideal gas at a specified temperature. When compared to average interparticle
spacing in the gas, it can be used to tell if the gas can be considered to be a classical or
Maxwell-Boltzmann gas, in which case the thermal wavelength must be much smaller than the average
interparticle spacing. Otherwise, quantum effects must be taken into account.

**Notation:**

#. :math:`\hbar` (:code:`hbar`) is the reduced Planck constant.
#. :math:`k_\text{B}` (:code:`k_B`) is the Boltzmann constant.
"""

from sympy import Eq, sqrt, pi
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
    clone_symbol,
)

thermal_wavelength = Symbol("thermal_wavelength", units.length)
r"""
Thermal de Broglies wavelength of the gas.

Symbol:
    :code:`lambda`

Latex:
    :math:`\lambda`
"""

particle_mass = clone_symbol(symbols.basic.mass, "particle_mass")
"""
:attr:`~symplyphysics.symbols.basic.mass` of a single gas particle.

Symbol:
    :code:`m`
"""

temperature = symbols.thermodynamics.temperature
"""
:attr:`~symplyphysics.symbols.thermodynamics.temperature` of the gas.

Symbol:
    :code:`T`
"""

definition = Eq(thermal_wavelength,
    units.hbar * sqrt(2 * pi / (particle_mass * units.boltzmann_constant * temperature)))
r"""
:code:`lambda = hbar * sqrt(2 * pi / (m * k_B * T))`

Latex:
    .. math::
        \lambda = \hbar \sqrt{\frac{2 \pi}{m k_\text{B} T}}
"""


@validate_input(
    particle_mass_=particle_mass,
    temperature_=temperature,
)
@validate_output(thermal_wavelength)
def calculate_thermal_wavelength(
    particle_mass_: Quantity,
    temperature_: Quantity,
) -> Quantity:
    result = definition.rhs.subs({
        particle_mass: particle_mass_,
        temperature: temperature_,
    })
    return Quantity(result)

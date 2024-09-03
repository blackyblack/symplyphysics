r"""
Spectral energy density at all frequencies
==========================================

*Planck's radiation law* describes the spectral density of electromagnetic radiation emitted
by a black body in thermal equlibrium at a given temperature when there is no net flow of
matter or energy between the body and its environment.

**Notation:**

#. :math:`h` is the Planck constant.
#. :math:`c` is the speed of light.
#. :math:`k_\text{B}` (:code:`k_B`) is the Boltzmann constant.

**Conditions:**

#. The black body is isolated from the environment.
"""

from sympy import Eq, exp, pi
from sympy.physics.units import planck, speed_of_light, boltzmann_constant
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
    clone_symbol,
)

spectral_energy_density = Symbol("spectral_energy_density",
    units.energy / (units.volume * units.frequency))
r"""
Spectral energy density, which is energy per unit volume per unit frequency.

Symbol:
    :code:`u_nu`

Latex:
    :math:`u_\nu`
"""

radiation_frequency = Symbol("radiation_frequency", units.frequency)
r"""
Frequency (linear) of the radiation.

Symbol:
    :code:`nu`

Latex:
    :math:`\nu`
"""

equilibrium_temperature = clone_symbol(symbols.thermodynamics.temperature)
"""
Equilibrium :attr:`~symplyphysics.symbols.thermodynamics.temperature` of the ensemble.
"""

law = Eq(spectral_energy_density, (8 * pi * planck * radiation_frequency**3 / speed_of_light**3) /
    (exp(planck * radiation_frequency / (boltzmann_constant * equilibrium_temperature)) - 1))
r"""
:code:`u_nu = (8 * pi * h * nu^3 / c^3) / (exp(h * nu / (k_B * T)) - 1)`

Latex:
    .. math::
        u_\nu = \frac{8 \pi h \nu^3}{c^3} \frac{1}{\exp \left( \frac{h \nu}{k_\text{B} T} \right) - 1}
"""


@validate_input(
    radiation_frequency_=radiation_frequency,
    equilibrium_temperature_=equilibrium_temperature,
)
@validate_output(spectral_energy_density)
def calculate_spectral_energy_density(
    radiation_frequency_: Quantity,
    equilibrium_temperature_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        radiation_frequency: radiation_frequency_,
        equilibrium_temperature: equilibrium_temperature_,
    })
    return Quantity(result)

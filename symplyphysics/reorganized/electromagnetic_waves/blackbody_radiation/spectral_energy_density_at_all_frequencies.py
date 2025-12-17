r"""
Spectral energy density at all frequencies
==========================================

*Planck's radiation law* describes the spectral density of electromagnetic radiation emitted
by a black body in thermal equilibrium at a given temperature when there is no net flow of
matter or energy between the body and its environment.

**Notation:**

#. :quantity_notation:`planck`.
#. :quantity_notation:`speed_of_light`.
#. :quantity_notation:`boltzmann_constant`.

**Conditions:**

#. The black body is isolated from the environment.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Planck%27s_law>`__.
"""

from sympy import Eq, exp, pi
from symplyphysics import Quantity, validate_input, validate_output, symbols
from symplyphysics.quantities import planck, speed_of_light, boltzmann_constant

spectral_energy_density = symbols.spectral_energy_density
"""
:symbols:`spectral_energy_density`.
"""

radiation_frequency = symbols.temporal_frequency
"""
:symbols:`temporal_frequency` of the radiation.
"""

equilibrium_temperature = symbols.temperature
"""
Equilibrium :symbols:`temperature` of the ensemble.
"""

law = Eq(spectral_energy_density, (8 * pi * planck * radiation_frequency**3 / speed_of_light**3) /
    (exp(planck * radiation_frequency / (boltzmann_constant * equilibrium_temperature)) - 1))
"""
:laws:symbol::

:laws:latex::
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

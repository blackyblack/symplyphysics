r"""
Spectral energy density at high frequency limit
===============================================

*Wien's approximation*, also known as *Wien distribution law*, describes the spectrum of blackbody
thermal radiation. It accurately describes short-wavelength (i.e. high-frequency) spectrum of
thermal emission, but fails to do that for long-wavelength (i.e. low-frequency) emission.

**Notation:**

#. :quantity_notation:`planck`.
#. :quantity_notation:`speed_of_light`.
#. :quantity_notation:`boltzmann_constant`.

**Conditions:**

#. The black body is isolated from the environment.
#. :math:`h \nu \gg k_\text{B} T`, i.e. photon energy is much greater than thermal energy.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Wien_approximation>`__.
"""

from sympy import Eq, exp, pi, Symbol as SymSymbol, solve, S
from symplyphysics import Quantity, validate_input, validate_output, symbols
from symplyphysics.quantities import planck, speed_of_light, boltzmann_constant
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.waves.blackbody_radiation import (
    spectral_energy_density_at_all_frequencies as planck_law,
)

spectral_energy_density = symbols.spectral_energy_density
"""
:symbols:`spectral_energy_density`.
"""

radiation_frequency = symbols.temporal_frequency
r"""
:symbols:`temporal_frequency` of the radiation.
"""

equilibrium_temperature = symbols.temperature
"""
Equilibrium :symbols:`temperature` of the ensemble.
"""

law = Eq(spectral_energy_density, (8 * pi * planck * radiation_frequency**3 / speed_of_light**3) *
    exp(-1 * planck * radiation_frequency / (boltzmann_constant * equilibrium_temperature)))
"""
:laws:symbol::

:laws:latex::
"""

# Derive from Planck's law of blackbody radiation

_planck_density = planck_law.law.rhs.subs({
    planck_law.radiation_frequency: radiation_frequency,
    planck_law.equilibrium_temperature: equilibrium_temperature,
})

_reduced_frequency = SymSymbol("reduced_frequency")

_reduced_frequency_eqn = Eq(
    _reduced_frequency,
    planck * radiation_frequency / (boltzmann_constant * equilibrium_temperature),
)

_planck_density_reduced = solve(
    (Eq(spectral_energy_density, _planck_density), _reduced_frequency_eqn),
    (spectral_energy_density, radiation_frequency),
    dict=True,
)[0][spectral_energy_density]

_planck_density_reduced_series = _planck_density_reduced.series(_reduced_frequency, S.Infinity,
    2).removeO()

_planck_density_series = solve(
    (Eq(spectral_energy_density, _planck_density_reduced_series), _reduced_frequency_eqn),
    (spectral_energy_density, _reduced_frequency),
    dict=True,
)[0][spectral_energy_density]

assert expr_equals(_planck_density_series, law.rhs)


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

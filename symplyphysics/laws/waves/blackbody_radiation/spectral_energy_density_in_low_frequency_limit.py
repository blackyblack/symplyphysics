r"""
Spectral energy density at low frequency limit
==============================================

The *Rayleigh-Jeans law* is an approximation to the spectral radiance of electromagnetic radiation
as a function of wave frequency from a blackbody at a given temperature through classical arguments.
The Rayleigh-Jeans law agrees with experimental results at large wavelengths (i.e. at low frequencies)
but strongly disagrees at short wavelengths (i.e. at high frequencies). This inconsistency is commonly
known as the *ultraviolet catastrophe*.

**Notation:**

#. :quantity_notation:`speed_of_light`.
#. :quantity_notation:`boltzmann_constant`.

**Conditions:**

#. The black body is isolated from the environment.
#. :math:`h \nu \ll k_\text{B} T`, i.e. photon energy is much smaller than thermal energy.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Rayleigh%E2%80%93Jeans_law>`__.
"""

from sympy import Eq, pi, Symbol as SymSymbol, solve
from symplyphysics import Quantity, validate_input, validate_output, symbols
from symplyphysics.quantities import planck, speed_of_light, boltzmann_constant
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.waves.blackbody_radiation import (
    spectral_energy_density_at_all_frequencies as planck_law,)

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

law = Eq(
    spectral_energy_density, 8 * pi * radiation_frequency**2 * boltzmann_constant *
    equilibrium_temperature / speed_of_light**3)
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

_planck_density_reduced_series = _planck_density_reduced.series(_reduced_frequency, 0, 3).removeO()

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

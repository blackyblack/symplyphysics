r"""
Spectral energy density at high frequency limit
===============================================

*Wien's approximation*, also known as *Wien distribution law*, describes the spectrum of blackbody
thermal radiation. It accurately describes short-wavelength (i.e. high-frequency) spectrum of
thermal emission, but fails to do that for long-wavelength (i.e. low-frequency) emission.

**Notation:**

#. :math:`h` is the Planck constant.
#. :math:`c` is the speed of light.
#. :math:`k_\text{B}` (:code:`k_B`) is the Boltzmann constant.

**Conditions:**

#. The black body is isolated from the environment.
#. :math:`h \nu \gg k_\text{B} T`, i.e. photon energy is much greater than thermal energy.
"""

from sympy import Eq, exp, pi, Symbol as SymSymbol, solve, S
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
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.waves.blackbody_radiation import spectral_energy_density_at_all_frequencies as planck_law

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

law = Eq(spectral_energy_density, (8 * pi * planck * radiation_frequency**3 / speed_of_light**3) *
    exp(-1 * planck * radiation_frequency / (boltzmann_constant * equilibrium_temperature)))
r"""
:code:`u_nu = (8 * pi * h * nu^3 / c^3) * exp(-1 * h * nu / (k_B * T))`

Latex:
    .. math::
        u_\nu = \frac{8 \pi h \nu^3}{c^3} \exp \left( - \frac{h \nu}{k_B T} \right)
"""

# Derive from Planck's law of blackbody radiation

_planck_density = planck_law.law.rhs.subs({
    planck_law.radiation_frequency: radiation_frequency,
    planck_law.equilibrium_temperature: equilibrium_temperature,
})

_reduced_frequency = SymSymbol("reduced_frequency")

_reduced_frequency_eqn = Eq(
    _reduced_frequency,
    units.planck * radiation_frequency / (boltzmann_constant * equilibrium_temperature),
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

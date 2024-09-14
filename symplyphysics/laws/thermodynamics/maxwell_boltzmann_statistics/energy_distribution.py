r"""
Energy distribution
===================

For a system containing a large number of identical non-interacting non-relativistic classical
particles in thermodynamic equilibrium, the energy distribution function is a function such that
:math:`f(E) dE` gives the fraction of particles with energies in the interval :math:`dE`
around energy value :math:`E`.

**Notation:**

#. :math:`k_\text{B}` (:code:`k_B`) is the Boltzmann constant.

**Notes:**

#. Number of particles is big enough that the laws of thermodynamics can be applied.
#. Particles are identical, non-interacting, non-relativistic, and classical.
#. The ensemble of particles is at thermodynamic equilibrium.
"""

from sympy import Eq, Rational, sqrt, pi, exp, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    clone_as_symbol,
    symbols,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.dynamics import kinetic_energy_from_mass_and_speed as kinetic_energy_law
from symplyphysics.laws.thermodynamics.maxwell_boltzmann_statistics import speed_distribution

energy_distribution_function = Symbol("energy_distribution_function",
    1 / units.energy,
    positive=True)
"""
Energy distribution function.

Symbol:
    :code:`f(E)`
"""

energy = Symbol("energy", units.energy, positive=True)
"""
Energy of the ensemble.

Symbol:
    :code:`E`
"""

equilibrium_temperature = clone_as_symbol(symbols.temperature, positive=True)
"""
Equilibrium :attr:`~symplyphysics.symbols.temperature` of the ensemble.
"""

law = Eq(
    energy_distribution_function,
    2 * sqrt(energy / pi) * (units.boltzmann_constant * equilibrium_temperature)**Rational(-3, 2) *
    exp(-1 * energy / (units.boltzmann_constant * equilibrium_temperature)))
r"""
:code:`f(E) = 2 * sqrt(E / pi) * (k_B * T)^(3/2) * exp(-1 * E / (k_B * T))`

Latex:
    .. math::
        f(E) = 2 \sqrt \frac{E}{\pi} \left( k_\text{B} T \right)^{3/2} \exp \left( - \frac{E}{k_\text{B} T} \right)
"""

# Derive from speed distribution and kinetic energy formula

# `E(v) = m*(v**2)/2` is a monotonous function of `v` for non-negative values of v.
# Therefore we can perform a substitution of variables in the distribution function
# via the [formula](https://en.wikipedia.org/wiki/Probability_density_function#Scalar_to_scalar)
# `f_E(E) = f_v(v(E)) * abs(dv(E)/dE)`

_speed = solve(kinetic_energy_law.law, kinetic_energy_law.speed)[0].subs({
    kinetic_energy_law.kinetic_energy: energy,
    kinetic_energy_law.mass: speed_distribution.particle_mass,
})

_speed_distribution = speed_distribution.law.rhs.subs(speed_distribution.equilibrium_temperature,
    equilibrium_temperature)

_speed_derivative_wrt_energy = _speed.diff(energy)

_energy_distribution_derived = (
    _speed_distribution.subs(speed_distribution.particle_speed, _speed) *
    abs(_speed_derivative_wrt_energy))

assert expr_equals(_energy_distribution_derived, law.rhs)


@validate_input(
    energy_=energy,
    equilibrium_temperature_=equilibrium_temperature,
)
@validate_output(energy_distribution_function)
def calculate_energy_distribution_function(
    energy_: Quantity,
    equilibrium_temperature_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        energy: energy_,
        equilibrium_temperature: equilibrium_temperature_,
    })
    return Quantity(result)

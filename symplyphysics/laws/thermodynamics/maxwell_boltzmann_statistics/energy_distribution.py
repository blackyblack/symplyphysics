from sympy import Eq, Rational, sqrt, pi, exp, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    Function,
    print_expression,
    validate_input,
    validate_output,
    clone_symbol,
    symbols,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.dynamics import kinetic_energy_from_mass_and_velocity as kinetic_energy_law
from symplyphysics.laws.thermodynamics.maxwell_boltzmann_statistics import speed_distribution

# For a system containing a large number of identical non-interacting non-relativistic classical
## particles in thermodynamic equilibrium, the energy distribution function is a function such that
## `f(E) * dE` gives the fraction of particles with energies in the interval `dE` around energy
## value `E`.

# Law: f(E) = 2 * sqrt(E / pi) * (k * T)**(-3/2) * exp(-E / (k * T))
## f(E) - energy distribution function
## E - energy
## k - Boltzmann constant
## T - equilibrium temperature

# Conditions
## - Number of particles is big enough that the laws of thermodynamics can be applied.
## - Particles are identical, non-interacting, non-relativistic, and obeying classical laws of physics.
## - The ensemble of particles is at thermodynamic equilibrium.

energy_distribution_function = Function("energy_distribution_function",
    1 / units.energy,
    positive=True)
energy = Symbol("energy", units.energy, positive=True)
equilibrium_temperature = clone_symbol(symbols.thermodynamics.temperature,
    "equilibrium_temperature",
    positive=True)

law = Eq(
    energy_distribution_function(energy),
    2 * sqrt(energy / pi) * (units.boltzmann_constant * equilibrium_temperature)**Rational(-3, 2) *
    exp(-1 * energy / (units.boltzmann_constant * equilibrium_temperature)))

# Derive from speed distribution and kinetic energy formula

# `E(v) = m*(v**2)/2` is a monotonous function of `v` for non-negative values of v.
# Therefore we can perform a substitution of variables in the distribution function
# via the [formula](https://en.wikipedia.org/wiki/Probability_density_function#Scalar_to_scalar)
# `f_E(E) = f_v(v(E)) * abs(dv(E)/dE)`

_speed = solve(kinetic_energy_law.law, kinetic_energy_law.body_velocity)[0].subs({
    kinetic_energy_law.kinetic_energy_of_body: energy,
    symbols.basic.mass: speed_distribution.particle_mass,
})

_speed_distribution = speed_distribution.law.rhs.subs(speed_distribution.equilibrium_temperature,
    equilibrium_temperature)

_speed_derivative_wrt_energy = _speed.diff(energy)

_energy_distribution_derived = (
    _speed_distribution.subs(speed_distribution.particle_speed, _speed) *
    abs(_speed_derivative_wrt_energy))

assert expr_equals(_energy_distribution_derived, law.rhs)


def print_law() -> str:
    return print_expression(law)


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

from sympy import Eq, Rational, sqrt, pi, exp
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

# For a system containing a large number of identical non-interacting non-relativistic classical
## particles in thermodynamic equilibrium, the speed distribution function is a function such that
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

energy_distribution_function = Function("energy_distribution_function", 1 / units.energy, positive=True)
energy = Symbol("energy", units.energy, positive=True)
equilibrium_temperature = clone_symbol(symbols.thermodynamics.temperature, "equilibrium_temperature", positive=True)

law = Eq(
    energy_distribution_function(energy),
    2
    * sqrt(energy / pi)
    * (units.boltzmann_constant * equilibrium_temperature)**Rational(-3, 2)
    * exp(-1 * energy / (units.boltzmann_constant * equilibrium_temperature))
)


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

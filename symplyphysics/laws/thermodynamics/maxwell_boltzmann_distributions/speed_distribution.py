from sympy import Eq, Rational, sqrt, pi, exp
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    clone_symbol,
    symbols,
)

# Description
## For a system containing a large number of identical non-interacting non-relativistic classical
## particles in thermodynamic equilibrium, the speed distribution function is a function such that
## `f(v) * dv` gives the fraction of particles with speeds in the interval `dv` at speed `v`.

# Law: f(v) = sqrt(2 / pi) * (m / (k * T))**(3/2) * v**2 * exp(-m * v**2 / (2 * k * T))
## f(v) - speed distribution function
## v - speed of particles, i.e. magnitude of velocity vector
## m - mass of particle
## k - Boltzmann constant
## T - equilibrium temperature of the particle ensemble

# Conditions
## - Number of particles is big enough that the laws of thermodynamics can be applied.
## - Particles are identical, non-interacting, non-relativistic, and obeying classical laws of physics.
## - The ensemble of particles is at thermodynamic equilibrium.

speed_distribution_function = Symbol("speed_distribution_function", 1 / units.velocity, positive=True)
particle_speed = Symbol("particle_speed", units.velocity, nonnegative=True)
particle_mass = clone_symbol(symbols.basic.mass, "particle_mass", positive=True)
equilibrium_temperature = clone_symbol(symbols.thermodynamics.temperature, "equilibrium_temperature", positive=True)

law = Eq(
    speed_distribution_function,
    sqrt(2 / pi)
    * (particle_mass / (units.boltzmann_constant * equilibrium_temperature))**Rational(3, 2)
    * particle_speed**2
    * exp(-1 * particle_mass * particle_speed**2 / (2 * units.boltzmann_constant * equilibrium_temperature))
)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    particle_speed_=particle_speed,
    particle_mass_=particle_mass,
    ensemble_temperature_=equilibrium_temperature,
)
@validate_output(speed_distribution_function)
def calculate_speed_distribution_function(
    particle_speed_: Quantity,
    particle_mass_: Quantity,
    ensemble_temperature_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        particle_speed: particle_speed_,
        particle_mass: particle_mass_,
        equilibrium_temperature: ensemble_temperature_,
    })
    return Quantity(result)

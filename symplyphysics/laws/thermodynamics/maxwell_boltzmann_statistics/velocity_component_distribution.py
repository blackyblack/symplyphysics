from sympy import Eq, sqrt, pi, exp
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

# Description
## For a system containing a large number of identical non-interacting non-relativistic
## classical particles in thermodynamic equilibrium, the velocity component distribution
## is a function `f(v_x)` such that `f(v_x) * dv_x` gives the fraction of particles with
## speeds in the interval `dv_x` around velocity component v_x. Also see Note below.

# Law: f(v_k) = sqrt(m / (2 * pi * k * T)) * exp(-m * v_k**2 / (2 * k * T))
## f(v_k) - distribution function of velocity component v_k
## v_k - component of velocity vector in Cartesian coordinates (k = x, y, z)
## m - mass of particle
## k - Boltzmann constant
## T - equilibrium temperature of the particle ensemble

# Note
## - Works for any velocity component in Cartesian coordinates: v_x, v_y, v_z.

# Conditions
## - Number of particles is big enough that the laws of thermodynamics can be applied.
## - Particles are identical, non-interacting, non-relativistic, and obeying classical laws of physics.
## - The ensemble of particles is at thermodynamic equilibrium.

velocity_component_distribution = Function("velocity_component_distribution", 1 / units.velocity, positive=True)
velocity_component = Symbol("velocity_component", units.velocity, positive=True)
particle_mass = clone_symbol(symbols.basic.mass, "particle_mass", positive=True)
equilibrium_temperature = clone_symbol(symbols.thermodynamics.temperature, "equilibrium_temperature", positive=True)

law = Eq(
    velocity_component_distribution(velocity_component),
    sqrt(particle_mass / (2 * pi * units.boltzmann_constant * equilibrium_temperature))
    * exp(-1 * particle_mass * velocity_component**2 / (2 * units.boltzmann_constant * equilibrium_temperature))
)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    velocity_component_=velocity_component,
    particle_mass_=particle_mass,
    ensemble_temperature_=equilibrium_temperature,
)
@validate_output(velocity_component_distribution)
def calculate_velocity_component_distribution(
    velocity_component_: Quantity,
    particle_mass_: Quantity,
    ensemble_temperature_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        velocity_component: velocity_component_,
        particle_mass: particle_mass_,
        equilibrium_temperature: ensemble_temperature_,
    })
    return Quantity(result)

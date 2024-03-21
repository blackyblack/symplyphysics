from sympy import Eq, sqrt, pi, integrate, S
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    symbols,
    clone_symbol,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics.maxwell_boltzmann_distributions import speed_distribution

# Description
## The average (mean) speed is the expected value of the speed distribution.

# Law: <v> = sqrt((8 / pi) * (k * T / m))
## <v> - average molecular speed
## k - Boltzmann constant
## T - equilibrium temperature
## m - molecular mass

# Conditions
## - Assuming the particles are distributed according to Maxwell-Boltzmann statistics.

average_speed = Symbol("average_speed", units.velocity, positive=True)
equilibrium_temperature = clone_symbol(symbols.thermodynamics.temperature, "equilibrium_temperature", positive=True)
particle_mass = clone_symbol(symbols.basic.mass, "particle_mass", positive=True)

law = Eq(
    average_speed, 
    sqrt(8 * units.boltzmann_constant * equilibrium_temperature / (pi * particle_mass)),
)

# Derive law from Maxwell-Boltzmann distribution function using the probability formula 
# of finding the average value of a function y(x) when the distribution function f(x) of 
# the random variable x is known, which is calculating the integral of y(x)*f(x) with the 
# integration limits being the lowest and highest possible values of x. In the case of 
# this law, it is the integral of v*f(v) with v ranging from 0 to infinity.

_distribution = speed_distribution.law.rhs.subs({
    speed_distribution.equilibrium_temperature: equilibrium_temperature,
    speed_distribution.particle_mass: particle_mass,
})

_average_speed_derived = integrate(
    speed_distribution.particle_speed * _distribution,
    (speed_distribution.particle_speed, 0, S.Infinity)
)

assert expr_equals(_average_speed_derived, law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    equilibrium_temperature_=equilibrium_temperature,
    particle_mass_=particle_mass,
)
@validate_output(average_speed)
def calculate_average_speed(
    equilibrium_temperature_: Quantity,
    particle_mass_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        equilibrium_temperature: equilibrium_temperature_,
        particle_mass: particle_mass_,
    })
    return Quantity(result)

from sympy import Eq, sqrt, pi, S, stats, Interval
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
from symplyphysics.laws.thermodynamics.maxwell_boltzmann_statistics import speed_distribution

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
equilibrium_temperature = clone_symbol(symbols.thermodynamics.temperature,
    "equilibrium_temperature",
    positive=True)
particle_mass = clone_symbol(symbols.basic.mass, "particle_mass", positive=True)

law = Eq(
    average_speed,
    sqrt(8 * units.boltzmann_constant * equilibrium_temperature / (pi * particle_mass)),
)

# Derive law from Maxwell-Boltzmann distribution function

_distribution = speed_distribution.law.rhs.subs({
    speed_distribution.equilibrium_temperature: equilibrium_temperature,
    speed_distribution.particle_mass: particle_mass,
})

_speed_random_variable = stats.ContinuousRV(
    speed_distribution.particle_speed,
    _distribution,
    set=Interval(0, S.Infinity),
)

# Average speed is the expected value of the speed random variable
_average_speed_derived = stats.E(_speed_random_variable)

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

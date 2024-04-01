from sympy import Eq, sqrt, solve, sign, S
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
## The most probable speed is the speed at which the Maxwell-Boltzmann speed distribution function
## is maximum.

# Law: v_prob = sqrt(2 * k * T / m)
## v_prob - most probable speed
## k - Boltzmann constant
## T - equilibrium temperature
## m - particle mass

most_probable_speed = Symbol("most_probable_speed", units.velocity, positive=True)
equilibrium_temperature = clone_symbol(symbols.thermodynamics.temperature,
    "equilibrium_temperature",
    positive=True)
particle_mass = clone_symbol(symbols.basic.mass, "particle_mass", positive=True)

law = Eq(most_probable_speed,
    sqrt(2 * units.boltzmann_constant * equilibrium_temperature / particle_mass))

# Derive from the Maxwell-Boltzmann speed distribution function

_distribution = speed_distribution.law.rhs.subs({
    speed_distribution.particle_mass: particle_mass,
    speed_distribution.equilibrium_temperature: equilibrium_temperature,
})

_distribution_first_derivative = _distribution.diff(speed_distribution.particle_speed)

# Found the points of extremum of the speed distribution function
_solutions = solve(_distribution_first_derivative, speed_distribution.particle_speed)

# Asserting there is only one point of extremum
assert len(_solutions) == 1
_most_probable_speed_derived = _solutions[0]

assert expr_equals(_most_probable_speed_derived, law.rhs)

_distribution_second_derivative = _distribution_first_derivative.diff(
    speed_distribution.particle_speed)

_distribution_second_derivative_at_most_probable_speed = _distribution_second_derivative.subs(
    speed_distribution.particle_speed, _most_probable_speed_derived)

# Proved that the point of extremum found is the point of maximum
assert sign(_distribution_second_derivative_at_most_probable_speed) == S.NegativeOne


def print_law() -> str:
    return print_expression(law)


@validate_input(
    equilibrium_temperature_=equilibrium_temperature,
    particle_mass_=particle_mass,
)
@validate_output(most_probable_speed)
def calculate_most_probable_speed(
    equilibrium_temperature_: Quantity,
    particle_mass_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        equilibrium_temperature: equilibrium_temperature_,
        particle_mass: particle_mass_,
    })
    return Quantity(result)

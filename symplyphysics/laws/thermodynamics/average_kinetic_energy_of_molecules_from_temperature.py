from sympy import Eq, solve, Rational, stats, Interval, S
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
from symplyphysics.laws.thermodynamics.maxwell_boltzmann_statistics import energy_distribution

# Description
## The kinetic theory of ideal gases allows us to determine the average kinetic energy for an ideal gas.
## It states that the average kinetic energy for all ideal gases is directly proportional to the absolute temperature of the gas and only depends on the temperature.

## Definition: E = 1.5 * k * T
## Where:
## E is average kinetic energy of molecules
## k is Boltzmann constant,
## T is temperature.

## Conditions
## The gas must be ideal

average_kinetic_energy = Symbol("average_kinetic_energy", units.energy, positive=True)
equilibrium_temperature = clone_symbol(symbols.thermodynamics.temperature, "equilibrium_temperature", positive=True)

law = Eq(average_kinetic_energy, Rational(3, 2) * units.boltzmann * equilibrium_temperature)

# Derive from Maxwell-Boltzmann energy distribution

_distribution = energy_distribution.law.rhs.subs(
    energy_distribution.equilibrium_temperature,
    equilibrium_temperature,
)

_random_energy_variable = stats.ContinuousRV(
    energy_distribution.energy,
    _distribution,
    set=Interval(0, S.Infinity),
)

_average_energy_derived = stats.E(_random_energy_variable)

assert expr_equals(_average_energy_derived, law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(temperature_=symbols.thermodynamics.temperature)
@validate_output(average_kinetic_energy)
def calculate_average_kinetic_energy(temperature_: Quantity) -> Quantity:
    result_expr = solve(law, average_kinetic_energy, dict=True)[0][average_kinetic_energy]
    result_average_kinetic_energy = result_expr.subs(equilibrium_temperature,
        temperature_)
    return Quantity(result_average_kinetic_energy)

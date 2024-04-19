from sympy import (Eq, solve, S, stats, Interval)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, symbols, clone_symbol)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics.maxwell_boltzmann_statistics import speed_distribution

# Description
## For an ideal gas, the average square of velocity is directly proportional to its temperature
## and inversely proportional to the molar mass of the gas.

# Law: <V^2> = 3 k T / m
## <V^2> - average square of velocity
## T - temperature in gas
## m - mass of molecule
## k - Boltzmann constant

average_square_velocity = Symbol("average_square_velocity", units.velocity**2, positive=True)
temperature_in_gas = clone_symbol(symbols.thermodynamics.temperature,
    "temperature_in_gas",
    positive=True)
molecule_mass = clone_symbol(symbols.basic.mass, "molecule_mass", positive=True)

law = Eq(average_square_velocity, 3 * units.boltzmann_constant * temperature_in_gas / molecule_mass)

# Derive law from Maxwell-Boltzmann distribution function

_distribution = speed_distribution.law.rhs.subs({
    speed_distribution.particle_mass: molecule_mass,
    speed_distribution.equilibrium_temperature: temperature_in_gas,
})

_speed_random_variable = stats.ContinuousRV(
    speed_distribution.particle_speed,
    _distribution,
    set=Interval(0, S.Infinity),
)

_average_square_of_speed_derived = stats.E(_speed_random_variable**2)

assert expr_equals(_average_square_of_speed_derived, law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(temperature_in_gas_=temperature_in_gas, mass_of_molecule_=molecule_mass)
@validate_output(average_square_velocity)
def calculate_average_square_velocity(temperature_in_gas_: Quantity,
    mass_of_molecule_: Quantity) -> Quantity:
    result_average_square_velocity = solve(law, average_square_velocity,
        dict=True)[0][average_square_velocity]
    result_expr = result_average_square_velocity.subs({
        temperature_in_gas: temperature_in_gas_,
        molecule_mass: mass_of_molecule_
    })
    return Quantity(result_expr)

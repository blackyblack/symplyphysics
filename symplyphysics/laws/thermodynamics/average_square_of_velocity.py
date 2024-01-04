from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression,
                           validate_input, validate_output)

# Description
## For an ideal gas, the average square of velocity is directly proportional to its temperature and inversely proportional to the molar mass of the gas: <V^2> = 3 k T / m
## Where:
## <V^2> - average square of velocity
## T - temperature in gas
## m - mass of a gas
## k - Boltzman's constant

temperature_in_gas = Symbol("temperature_in_gas", units.temperature)
mass_of_gas = Symbol("mass_of_gas", units.mass)

average_square_velocity = Symbol("average_square_velocity", units.velocity ** 2)

law = Eq(average_square_velocity, 3 * units.boltzmann_constant * temperature_in_gas / mass_of_gas)


def print_law() -> str:
    return print_expression(law)


@validate_input(temperature_in_gas_=temperature_in_gas, molar_mass_=mass_of_gas)
@validate_output(average_square_velocity)
def calculate_average_square_velocity(temperature_in_gas_: Quantity, mass_of_gas_: Quantity) -> Quantity:
    result_average_square_velocity = solve(law, average_square_velocity, dict=True)[0][average_square_velocity]
    result_expr = result_average_square_velocity.subs({temperature_in_gas: temperature_in_gas_, mass_of_gas: mass_of_gas_})
    return Quantity(result_expr)

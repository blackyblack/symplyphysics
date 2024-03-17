from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, symbols, clone_symbol)

# Description
## For an ideal gas, the average square of velocity is directly proportional to its temperature and inversely proportional to the molar mass of the gas: <V^2> = 3 k T / m
## Where:
## <V^2> - average square of velocity
## T - temperature in gas
## m - mass of molecule
## k - Boltzman's constant

temperature_in_gas = Symbol("temperature_in_gas", units.temperature)

average_square_velocity = Symbol("average_square_velocity", units.velocity**2)
molecule_mass = clone_symbol(symbols.basic.mass, "molecule_mass")

law = Eq(average_square_velocity, 3 * units.boltzmann_constant * temperature_in_gas / molecule_mass)


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

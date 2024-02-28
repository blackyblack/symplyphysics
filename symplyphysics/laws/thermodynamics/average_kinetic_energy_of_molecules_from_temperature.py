from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, symbols)

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

average_kinetic_energy = Symbol("average_kinetic_energy", units.energy)

law = Eq(average_kinetic_energy, 1.5 * units.boltzmann * symbols.basic.temperature)


def print_law() -> str:
    return print_expression(law)


@validate_input(temperature_=symbols.basic.temperature)
@validate_output(average_kinetic_energy)
def calculate_average_kinetic_energy(temperature_: Quantity) -> Quantity:
    result_expr = solve(law, average_kinetic_energy, dict=True)[0][average_kinetic_energy]
    result_average_kinetic_energy = result_expr.subs(symbols.basic.temperature, temperature_)
    return Quantity(result_average_kinetic_energy)

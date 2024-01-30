from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)

# Description

## Ideal gas equation: p = 2/3 * n * E
## Where:
## p is pressure
## n is concentration of molecules
## E is the average kinetic energy of gas molecules

## Conditions
## The gas must be ideal

average_kinetic_energy = Symbol("average_kinetic_energy", units.energy)
molecules_concentration = Symbol("molecules_concentration", 1 / units.volume)
pressure = Symbol("pressure", units.pressure)

law = Eq(pressure, 2 / 3 * molecules_concentration * average_kinetic_energy)

def print_law() -> str:
    return print_expression(law)


@validate_input(molecules_concentration_=molecules_concentration, average_kinetic_energy_=average_kinetic_energy )
@validate_output(pressure)
def calculate_pressure(molecules_concentration_: Quantity, average_kinetic_energy_: Quantity) -> Quantity:
    result_expr = solve(law, pressure, dict=True)[0][pressure]
    result_pressure = result_expr.subs({
        molecules_concentration: molecules_concentration_,
        average_kinetic_energy: average_kinetic_energy_,
    })
    return Quantity(result_pressure)

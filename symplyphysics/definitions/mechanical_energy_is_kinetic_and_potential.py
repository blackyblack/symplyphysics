from sympy import (Eq, solve)
from symplyphysics import (units, expr_to_quantity, Quantity, Symbol, print_expression,
    validate_input_symbols, validate_output_symbol)

# Description
## The mechanical energy of the system is defined as the total kinetic energy plus the total potential energy.

# Definition: E = K + P
# Where:
## K is kinetic energy of a system,
## P is potential energy of a system,
## E is mechanical energy.

mechanical_energy = Symbol("mechanical_energy", units.energy)
kinetic_energy = Symbol("kinetic_energy", units.energy)
potential_energy = Symbol("potential_energy", units.energy)

definition = Eq(mechanical_energy, kinetic_energy + potential_energy)

definition_units_SI = units.joule


def print() -> str:
    return print_expression(definition)


@validate_input_symbols(kinetic_energy_=kinetic_energy, potential_energy_=potential_energy)
@validate_output_symbol(mechanical_energy)
def calculate_mechanical_energy(kinetic_energy_: Quantity, potential_energy_: Quantity) -> Quantity:
    solved = solve(definition, mechanical_energy, dict=True)[0][mechanical_energy]
    result_expr = solved.subs({
        kinetic_energy: kinetic_energy_,
        potential_energy: potential_energy_
    })
    return expr_to_quantity(result_expr)

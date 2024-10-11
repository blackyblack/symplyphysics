from sympy import Eq, solve
from symplyphysics import (Quantity, Symbol, validate_input, validate_output, dimensionless,
    convert_to_float, clone_as_symbol, symbols)

# Description
## The traveling atom moves towards the substrate in the magnetron. At the same time, it collides with gas atoms.
## The energy transfer coefficient in these collisions depends on the mass of the traveling atom and the mass of the gas atom.

## Law is: x = 2 * M1 * M2 / (M1 + M2)^2, where
## x - energy transfer coefficient,
## M1 - the mass of the traveling atom,
## M2 - the mass of the gas atom.

energy_transfer_coefficient = Symbol("energy_transfer_coefficient", dimensionless)

mass_of_traveling_atom = clone_as_symbol(symbols.mass, display_symbol="m_1", display_latex="m_1")
mass_of_gas_atom = clone_as_symbol(symbols.mass, display_symbol="m_2", display_latex="m_2")

law = Eq(
    energy_transfer_coefficient,
    2 * mass_of_traveling_atom * mass_of_gas_atom / (mass_of_traveling_atom + mass_of_gas_atom)**2)


@validate_input(mass_of_traveling_atom_=mass_of_traveling_atom, mass_of_gas_atom_=mass_of_gas_atom)
@validate_output(energy_transfer_coefficient)
def calculate_energy_transfer_coefficient(mass_of_traveling_atom_: Quantity,
    mass_of_gas_atom_: Quantity) -> float:
    result_expr = solve(law, energy_transfer_coefficient, dict=True)[0][energy_transfer_coefficient]
    result_expr = result_expr.subs({
        mass_of_traveling_atom: mass_of_traveling_atom_,
        mass_of_gas_atom: mass_of_gas_atom_,
    })
    return convert_to_float(result_expr)

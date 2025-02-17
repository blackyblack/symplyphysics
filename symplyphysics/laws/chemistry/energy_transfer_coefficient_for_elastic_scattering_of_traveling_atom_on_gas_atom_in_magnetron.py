"""
Energy transfer coefficient for elastic scattering in magnetron
===============================================================

The traveling atom moves towards the substrate in the magnetron. At the same time, it
collides with gas atoms. The energy transfer coefficient in these collisions depends on
the mass of the traveling atom and the mass of the gas atom.

..
    TODO: find link
"""

from sympy import Eq, solve
from symplyphysics import (Quantity, Symbol, validate_input, validate_output, dimensionless,
    convert_to_float, clone_as_symbol, symbols)

energy_transfer_coefficient = Symbol("x", dimensionless)
"""
Energy transfer coefficient.
"""

traveling_atom_mass = clone_as_symbol(symbols.mass, subscript="1")
"""
:symbols:`mass` of traveling atom.
"""

gas_atom_mass = clone_as_symbol(symbols.mass, subscript="2")
"""
:symbols:`mass` of gas atom.
"""

law = Eq(energy_transfer_coefficient,
    2 * traveling_atom_mass * gas_atom_mass / (traveling_atom_mass + gas_atom_mass)**2)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(mass_of_traveling_atom_=traveling_atom_mass, mass_of_gas_atom_=gas_atom_mass)
@validate_output(energy_transfer_coefficient)
def calculate_energy_transfer_coefficient(mass_of_traveling_atom_: Quantity,
    mass_of_gas_atom_: Quantity) -> float:
    result_expr = solve(law, energy_transfer_coefficient, dict=True)[0][energy_transfer_coefficient]
    result_expr = result_expr.subs({
        traveling_atom_mass: mass_of_traveling_atom_,
        gas_atom_mass: mass_of_gas_atom_,
    })
    return convert_to_float(result_expr)

"""
Electrostatic potential energy of two charges via distance
==========================================================

Electrostatic potential energy due to two point charges depends on the inverse 
distance to the distance between the charges. Note that this is the energy of
interaction belonging to the entire system.
"""

from sympy import (Eq, solve, pi)
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

electrostatic_potential_energy = clone_as_symbol(symbols.energy, display_symbol="U_E", display_latex="U_\\mathbf{E}")
"""
Electrostatic potential :symbols:`energy` of system.
"""

absolute_permittivity = symbols.absolute_permittivity
"""
:symbols:`absolute_permittivity` of the medium.
"""

distance = symbols.distance_to_origin
"""
:symbols:`distance` between the point charges.
"""

first_charge = clone_as_symbol(symbols.charge, display_symbol="q_1", display_latex="q_1")
"""
Value of the first :symbols:`charge`.
"""

second_charge = clone_as_symbol(symbols.charge, display_symbol="q_2", display_latex="q_2")
"""
Value of the second :symbols:`charge`.
"""

law = Eq(electrostatic_potential_energy,
    (first_charge * second_charge) / (4 * pi * absolute_permittivity * distance))
"""
:laws:symbol::

:laws:latex::
"""

@validate_input(absolute_permittivity_=absolute_permittivity,
    distance_=distance,
    charge_1_=first_charge,
    charge_2_=second_charge)
@validate_output(electrostatic_potential_energy)
def calculate_energy(absolute_permittivity_: Quantity, distance_: Quantity, charge_1_: Quantity,
    charge_2_: Quantity) -> Quantity:
    result_expr = solve(law, electrostatic_potential_energy, dict=True)[0][electrostatic_potential_energy]
    result_expr = result_expr.subs({
        absolute_permittivity: absolute_permittivity_,
        distance: distance_,
        first_charge: charge_1_,
        second_charge: charge_2_,
    })
    return Quantity(result_expr)

"""
Magnetic field of coil
======================

Near the center of the coil, the magnetic field is quite uniform and directly
proportional to the current in the coil's wire.

**Notation:**

#. :quantity_notation:`vacuum_permeability`.

**Conditions:**

#. The medium is vacuum.
#. The magnetic field is measured near the center of the coil.

**Links:**

#. `Physics LibreTexts, similar formula 12.7.10 <https://phys.libretexts.org/Bookshelves/University_Physics/University_Physics_(OpenStax)/University_Physics_II_-_Thermodynamics_Electricity_and_Magnetism_(OpenStax)/12%3A_Sources_of_Magnetic_Fields/12.07%3A_Solenoids_and_Toroids>`__.
"""

from sympy import Eq, solve
from symplyphysics import Quantity, validate_input, validate_output, symbols, quantities

magnetic_flux_density = symbols.magnetic_flux_density
"""
:symbols:`magnetic_flux_density`.
"""

current = symbols.current
"""
:symbols:`current` flowing through the coil.
"""

length = symbols.length
"""
:symbols:`length` of the coil.
"""

coil_turn_count = symbols.positive_number
"""
Number of turns in the coil. See :symbols:`positive_number`.
"""

law = Eq(magnetic_flux_density, quantities.vacuum_permeability * current * coil_turn_count / length)
"""
:laws:symbol::

:laws:latex::
"""

# Refer to the link for the derivation.


@validate_input(current_=current, length_=length, number_turns_=coil_turn_count)
@validate_output(magnetic_flux_density)
def calculate_induction(current_: Quantity, length_: Quantity, number_turns_: float) -> Quantity:
    if number_turns_ < 0:
        raise ValueError("Number of turns cannot be negative")
    result_expr = solve(law, magnetic_flux_density, dict=True)[0][magnetic_flux_density]
    result_expr = result_expr.subs({
        current: current_,
        length: length_,
        coil_turn_count: number_turns_,
    })
    return Quantity(result_expr)

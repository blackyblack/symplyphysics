"""
Inductance via number of turns and coil volume
==============================================

The inductance of a coil (a solenoid) can be expressed as a function of the material's
permeability, the number of turns in the coil per unit length and its volume.

**Links:**

#. `Physics LibreTexts, formula 14.3.12 <https://phys.libretexts.org/Bookshelves/University_Physics/University_Physics_(OpenStax)/University_Physics_II_-_Thermodynamics_Electricity_and_Magnetism_(OpenStax)/14%3A_Inductance/14.03%3A_Self-Inductance_and_Inductors>`__.

..
    TODO rename file
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    SymbolNew,
    units,
)

inductance = symbols.inductance
"""
:symbols:`inductance` of the coil.
"""

absolute_permeability = symbols.absolute_permeability
"""
:symbols:`absolute_permeability` of the inside of the coil.
"""

specific_coil_turn_count = SymbolNew("n", 1 / units.length)
"""
Number of turns in the coil per unit length.
"""

volume = symbols.volume
"""
:symbols:`volume` of the coil.
"""

law = Eq(inductance,
    absolute_permeability * specific_coil_turn_count**2 * volume)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(absolute_permeability_=absolute_permeability,
    turn_count_=specific_coil_turn_count,
    volume_=volume)
@validate_output(inductance)
def calculate_inductance(absolute_permeability_: Quantity, turn_count_: Quantity,
    volume_: Quantity) -> Quantity:
    result_inductance_expr = solve(law, inductance, dict=True)[0][inductance]
    result_expr = result_inductance_expr.subs({
        absolute_permeability: absolute_permeability_,
        specific_coil_turn_count: turn_count_,
        volume: volume_
    })
    return Quantity(result_expr)

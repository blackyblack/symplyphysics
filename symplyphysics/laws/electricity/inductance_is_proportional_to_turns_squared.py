"""
Inductance is proportional to turn count
========================================

The basic characteristic of a coil is its inductance, which is the ability of the coil
to accumulate energy as magnetic field.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Inductance#Inductance_of_a_solenoid>`__.

..
    TODO: replace `mu_0 * mu_r` with `mu`
    TODO: fix file name
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
    dimensionless,
    quantities,
    symbols,
)

inductance = symbols.inductance
"""
:symbols:`inductance` of the coil.
"""

relative_permeability = symbols.relative_permeability
"""
:symbols:`relative_permeability` of the medium within the coil.
"""

turn_count = SymbolNew("N", dimensionless)
"""
Number of turns in the coil.
"""

cross_sectional_area = symbols.area
"""
Cross sectional :symbols:`area` of the coil.
"""

length = symbols.length
"""
:symbols:`length` of the coil.
"""

law = Eq(inductance,
    quantities.vacuum_permeability * relative_permeability * turn_count**2 * cross_sectional_area / length)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(turn_area_=cross_sectional_area, coil_length_=length)
@validate_output(inductance)
def calculate_inductance(magnetic_permeability_: float, number_of_turns_: float,
    turn_area_: Quantity, coil_length_: Quantity) -> Quantity:
    result_inductance_expr = solve(law, inductance, dict=True)[0][inductance]
    result_expr = result_inductance_expr.subs({
        relative_permeability: magnetic_permeability_,
        turn_count: number_of_turns_,
        cross_sectional_area: turn_area_,
        length: coil_length_
    })
    return Quantity(result_expr)

"""
Number of impacts on the wall from area and speed
=================================================

The number of impacts of molecules against the wall in a specified direction is
proportional to the number density of the molecules, the area of the wall, the
projection of velocity in said direction and time.

**Conditions:**

#. Gas is ideal.
#. Wall is flat and perpendicular to specified axis.

**Links:**

#. `Chemistry LibreTexts, similar formula <https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Physical_Chemistry_(LibreTexts)/27%3A_The_Kinetic_Theory_of_Gases/27.04%3A_The_Frequency_of_Collisions_with_a_Wall>`__.

..
    TODO: fix file name
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
    dimensionless,
    convert_to_float,
    symbols,
)

impact_count = SymbolNew("N", dimensionless)
"""
Impact count.
"""

number_density = symbols.number_density
"""
:symbols:`number_density` of particles.
"""

area = symbols.area
"""
Wall :symbols:`area`.
"""

velocity_projection = symbols.speed
"""
Projection of velocity vector in the given direction. See :symbols:`speed`.
"""

time = symbols.time
"""
:symbols:`time`.
"""

law = Eq(impact_count, (number_density * area * velocity_projection * time) / 2)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(molecules_concentration_=number_density,
    area_=area,
    velocity_projection_=velocity_projection,
    time_=time)
@validate_output(impact_count)
def calculate_number_of_impacts(molecules_concentration_: Quantity, area_: Quantity,
    velocity_projection_: Quantity, time_: Quantity) -> int:
    result_expr = solve(law, impact_count, dict=True)[0][impact_count]
    result_number_of_impacts = result_expr.subs({
        number_density: molecules_concentration_,
        area: area_,
        velocity_projection: velocity_projection_,
        time: time_
    })
    return int(convert_to_float(result_number_of_impacts))

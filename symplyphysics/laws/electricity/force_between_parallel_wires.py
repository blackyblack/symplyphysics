"""
Force between parallel wires
============================

Two parallel wires through which current flows interact with each other. The intensity
of their interaction depends on the magnitude of the current, the distance from the
wires, as well as on the material and length of the wire.

**Notes:**

#. Also see `Ampère's force law <https://en.wikipedia.org/wiki/Amp%C3%A8re%27s_force_law>`__.

**Links:**

#. `Physics LibreTexts, formula 22.10.3 <https://phys.libretexts.org/Bookshelves/College_Physics/College_Physics_1e_(OpenStax)/22%3A_Magnetism/22.10%3A_Magnetic_Force_between_Two_Parallel_Conductors>`__.
"""

from sympy import Eq, solve, pi
from symplyphysics import (
    symbols,
    Quantity,
    validate_input,
    validate_output,
    clone_as_symbol,
)

force = symbols.force
"""
:symbols:`force` of interaction between the wires.
"""

absolute_permeability = symbols.absolute_permeability
"""
:symbols:`absolute_permeability` of the medium.
"""

first_current = clone_as_symbol(symbols.current, subscript="1")
"""
:symbols:`current` in the first wire.
"""

second_current = clone_as_symbol(symbols.current, subscript="2")
"""
:symbols:`current` in the second wire.
"""

length = symbols.length
"""
:symbols:`length` of the wires.
"""

distance = symbols.euclidean_distance
"""
Perpendicular :symbols:`euclidean_distance` between the wires.
"""

law = Eq(
    force, absolute_permeability * first_current *
    second_current * length / (2 * pi * distance))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(absolute_permeability_=absolute_permeability,
    first_wire_current_=first_current,
    second_wire_current_=second_current,
    length_=length,
    distance_=distance)
@validate_output(force)
def calculate_force(absolute_permeability_: Quantity, first_wire_current_: Quantity,
    second_wire_current_: Quantity, length_: Quantity, distance_: Quantity) -> Quantity:
    result_expr = solve(law, force, dict=True)[0][force]
    result_expr = result_expr.subs({
        absolute_permeability: absolute_permeability_,
        first_current: first_wire_current_,
        second_current: second_wire_current_,
        length: length_,
        distance: distance_
    })
    return Quantity(result_expr)

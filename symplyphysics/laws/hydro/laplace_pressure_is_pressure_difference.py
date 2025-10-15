"""
Laplace pressure is pressure difference
=======================================

The **Laplace pressure** is the pressure difference between the inside and the outside of a curved
surface that forms the boundary between two fluid regions. This is caused by the surface tension at
the fluid interface.

The terms *outside* and *inside* are related to the direction dictated by the curvature of the
interface, not by the position relative to the fluid in question.

**Links:**

#. `Wikipedia — Laplace pressure <https://en.wikipedia.org/wiki/Laplace_pressure>`__.
"""

from sympy import Eq
from symplyphysics import Quantity, validate_input, validate_output, symbols, clone_as_symbol

laplace_pressure = clone_as_symbol(
    symbols.pressure,
    display_symbol="p_L",
    display_latex="p_\\text{L}",
)
"""
Laplace :symbols:`pressure` of the fluid boundary.
"""

outside_pressure = clone_as_symbol(
    symbols.pressure,
    display_symbol="p_o",
    display_latex="p_\\text{o}",
)
"""
:symbols:`pressure` at the side that the interface is bulging *from*.
"""

inside_pressure = clone_as_symbol(
    symbols.pressure,
    display_symbol="p_i",
    display_latex="p_\\text{i}",
)
"""
:symbols:`pressure` at the side that the interface is bulging *toward*.
"""

law = Eq(laplace_pressure, outside_pressure - inside_pressure)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    external_pressure_=outside_pressure,
    internal_pressure_=inside_pressure,
)
@validate_output(laplace_pressure)
def calculate_laplace_pressure(
    external_pressure_: Quantity,
    internal_pressure_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        outside_pressure: external_pressure_,
        inside_pressure: internal_pressure_,
    })

    return Quantity(result)


# UNIQUE_LAW_ID: 734

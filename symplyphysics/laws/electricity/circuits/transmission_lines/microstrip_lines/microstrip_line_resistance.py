"""
Resistance of microstrip line
=============================

The microstrip line is a dielectric substrate on which a metal strip is applied. Its
resistance depends on the physical dimensions of the microstrip and the resistance of
its surface.

..
    TODO: find link
"""

from sympy import Eq, solve, log
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

resistance = symbols.electrical_resistance
"""
:symbols:`electrical_resistance` of the microstrip line.
"""

length = symbols.length
"""
:symbols:`length` of the microstrip.
"""

width = clone_as_symbol(symbols.length, display_symbol="w", display_latex="w")
"""
Width (see :symbols:`length`) of the microstrip.
"""

thickness = clone_as_symbol(symbols.thickness, display_symbol="t", display_latex="t")
"""
:symbols:`thickness` of the microstrip.
"""

surface_resistance = clone_as_symbol(symbols.electrical_resistance, display_symbol="R_s", display_latex="R_\\text{s}")
"""
:symbols:`electrical_resistance` of the surface of the microstrip.
"""

law = Eq(
    resistance,
    (1.4 + 0.217 * log(width / (5 * thickness)))
    * (surface_resistance * length / (2 * (width + thickness))),
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(strip_thickness_=thickness,
    strip_length_=length,
    strip_width_=width,
    surface_resistance_=surface_resistance)
@validate_output(resistance)
def calculate_resistance(strip_thickness_: Quantity, strip_length_: Quantity,
    strip_width_: Quantity, surface_resistance_: Quantity) -> Quantity:
    result_expr = solve(law, resistance, dict=True)[0][resistance]
    result_expr = result_expr.subs({
        thickness: strip_thickness_,
        length: strip_length_,
        width: strip_width_,
        surface_resistance: surface_resistance_,
    })
    return Quantity(result_expr)

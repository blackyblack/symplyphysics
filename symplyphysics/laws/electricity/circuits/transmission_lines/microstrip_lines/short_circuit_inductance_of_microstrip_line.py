"""
Short circuit inductance of microstrip line
===========================================

The microstrip line is a dielectric substrate on which a metal strip is applied. To
create a short circuit based on a microstrip line, a metallized hole can be made in the
line. Such an opening will have a certain inductance, which can be calculated.

**Notation:**

#. :quantity_notation:`vacuum_permeability`.

..
    TODO: add link
"""

from sympy import Eq, solve, sqrt, log, pi, evaluate
from symplyphysics import Quantity, validate_input, validate_output, quantities, symbols

inductance = symbols.inductance
"""
:symbols:`inductance` of a metallized hole in a microstrip line.
"""

substrate_thickness = symbols.thickness
"""
:symbols:`thickness` of the substrate of the microstrip line.
"""

radius = symbols.radius
"""
:symbols:`radius` of the metallized hole.
"""

# the following block prevents the re-ordering of terms for the code printer
with evaluate(False):
    _constant = quantities.vacuum_permeability / (2 * pi)
    _first_expression = sqrt(radius**2 + substrate_thickness**2)
    _second_expression = substrate_thickness * log((substrate_thickness + _first_expression) / radius)
    _third_expression = 1.5 * (radius - _first_expression)

law = Eq(inductance, _constant * (_second_expression + _third_expression))
"""
:laws:symbol::

:laws:latex::

..
    TODO: check if `1.5` isn't really `3/2`
"""


@validate_input(
    thickness_of_substrate_=substrate_thickness,
    radius_of_hole_=radius,
)
@validate_output(inductance)
def calculate_inductance(thickness_of_substrate_: Quantity, radius_of_hole_: Quantity) -> Quantity:
    result_expr = solve(law, inductance, dict=True)[0][inductance]
    result_expr = result_expr.subs({
        substrate_thickness: thickness_of_substrate_,
        radius: radius_of_hole_,
    })
    return Quantity(result_expr)

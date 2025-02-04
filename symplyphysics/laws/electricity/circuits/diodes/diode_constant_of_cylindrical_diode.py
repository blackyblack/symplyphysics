"""
Diode constant of cylindrical diode
===================================

In a cylindrical diode, the cathode is located in the center, and the anode is located
around it in the form of a cylinder. The radius of the anode is usually much larger than
the radius of the cathode. The diode constant depends on the radii and on the area of
the anode.

**Notation:**

#. :quantity_notation:`vacuum_permittivity`.
#. :quantity_notation:`elementary_charge`.
#. :quantity_notation:`electron_rest_mass`.

..
    TODO: find link
"""

from sympy import Eq, Rational, solve, sqrt
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)
from symplyphysics.quantities import vacuum_permittivity, elementary_charge, electron_rest_mass

diode_constant = symbols.diode_constant
"""
:symbols:`diode_constant`.
"""

anode_area = clone_as_symbol(symbols.area, display_symbol="A_a", display_latex="A_\\text{a}")
"""
:symbols:`area` of the anode.
"""

anode_radius = clone_as_symbol(symbols.radius, display_symbol="r_a", display_latex="r_\\text{a}")
"""
:symbols:`radius` of the anode.
"""

cathode_radius = clone_as_symbol(symbols.radius, display_symbol="r_c", display_latex="r_\\text{c}")
"""
:symbols:`radius` of the cathode.
"""

law = Eq(
    diode_constant,
    Rational(4, 9) * vacuum_permittivity * sqrt(2 * elementary_charge / electron_rest_mass) *
    anode_area / (anode_radius**2 * (1 - cathode_radius / anode_radius)**2))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(anode_area_=anode_area, anode_radius_=anode_radius, cathode_radius_=cathode_radius)
@validate_output(diode_constant)
def calculate_diode_constant(anode_area_: Quantity, anode_radius_: Quantity,
    cathode_radius_: Quantity) -> Quantity:
    if anode_radius_.scale_factor <= cathode_radius_.scale_factor:
        raise ValueError("The anode radius must be greater than the cathode radius")
    result_expr = solve(law, diode_constant, dict=True)[0][diode_constant]
    result_expr = result_expr.subs({
        anode_area: anode_area_,
        anode_radius: anode_radius_,
        cathode_radius: cathode_radius_,
    })
    return Quantity(result_expr)

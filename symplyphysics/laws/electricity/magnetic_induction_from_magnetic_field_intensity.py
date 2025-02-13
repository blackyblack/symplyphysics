"""
Magnetic flux density from magnetic field strength
==================================================

Magnetic flux density is a physical quantity that is a force characteristic of a magnetic
field, namely, a characteristic of its action on moving charged particles and on bodies
with a magnetic moment.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Magnetic_field#The_H-field>`__.

..
    TODO: rename file
"""

from sympy import Eq, solve
from symplyphysics import Quantity, validate_input, validate_output, symbols

magnetic_flux_density = symbols.magnetic_flux_density
"""
:symbols:`magnetic_flux_density`.
"""

absolute_permeability = symbols.absolute_permeability
"""
:symbols:`absolute_permeability` of the medium.
"""

magnetic_field_strength = symbols.magnetic_field_strength
"""
:symbols:`magnetic_field_strength`.
"""

law = Eq(magnetic_flux_density, absolute_permeability * magnetic_field_strength)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(absolute_permeability_=absolute_permeability, intensity_=magnetic_field_strength)
@validate_output(magnetic_flux_density)
def calculate_induction(absolute_permeability_: Quantity, intensity_: Quantity) -> Quantity:
    result_expr = solve(law, magnetic_flux_density, dict=True)[0][magnetic_flux_density]
    result_expr = result_expr.subs({
        absolute_permeability: absolute_permeability_,
        magnetic_field_strength: intensity_,
    })
    return Quantity(result_expr)

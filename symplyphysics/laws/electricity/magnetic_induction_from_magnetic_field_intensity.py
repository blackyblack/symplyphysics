"""
Magnetic flux density from magnetic field strength
==================================================

Magnetic flux density is a physical quantity that is a force characteristic of a magnetic
field, namely, a characteristic of its action on moving charged particles and on bodies
with a magnetic moment.

**Notation:**

#. :quantity_notation:`vacuum_permeability`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Magnetic_field#The_H-field>`__.

..
    TODO: rename file
    TODO: replace `mu_0 * mu_r` with `mu`
"""

from sympy import Eq, solve
from symplyphysics import Quantity, validate_input, validate_output, quantities, symbols

magnetic_flux_density = symbols.magnetic_flux_density
"""
:symbols:`magnetic_flux_density`.
"""

relative_permeability = symbols.relative_permeability
"""
:symbols:`relative_permeability` of the medium.
"""

magnetic_field_strength = symbols.magnetic_field_strength
"""
:symbols:`magnetic_field_strength`.
"""

law = Eq(magnetic_flux_density, quantities.vacuum_permeability * relative_permeability * magnetic_field_strength)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(relative_permeability_=relative_permeability, intensity_=magnetic_field_strength)
@validate_output(magnetic_flux_density)
def calculate_induction(relative_permeability_: float, intensity_: Quantity) -> Quantity:
    result_expr = solve(law, magnetic_flux_density, dict=True)[0][magnetic_flux_density]
    result_expr = result_expr.subs({
        relative_permeability: relative_permeability_,
        magnetic_field_strength: intensity_,
    })
    return Quantity(result_expr)

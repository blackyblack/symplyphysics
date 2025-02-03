r"""
Submerged volume of floating body via density ratio
===================================================

If a body is fully or partially submerged in a fluid, an Archimedes force
starts acting on it from the fluid pushing it in the opposite direction of
the gravity force.

The ratio of the body's volume submerged in a fluid to its total volume
depends on the ratio of the densities of the body and the fluid.

**Conditions:**

#. :math:`\rho \le \rho_\text{fl}`, so the body must be floating. See below for the
   description of the symbols.
#. The body must be in static equilibrium.

**Links:**

#. `Physics LibreTexts, formula 10.3.17 <https://phys.libretexts.org/Bookshelves/University_Physics/Physics_(Boundless)/10%3A_Fluids/10.3%3A_Archimedes_Principle>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    convert_to_si,
    symbols,
    clone_as_symbol,
)

submerged_volume = clone_as_symbol(symbols.volume, display_symbol="V_fl", display_latex="V_\\text{fl}")
"""
:symbols:`volume` submerged in the fluid, which is equal to the volume of the displaced
fluid.
"""

body_volume = symbols.volume
"""
Total :symbols:`volume` of the body.
"""

body_density = symbols.density
"""
:symbols:`density` of the body.
"""

fluid_density = clone_as_symbol(symbols.density, display_symbol="rho_fl", display_latex="\\rho_\\text{fl}")
"""
:symbols:`density` of the fluid.
"""

law = Eq(submerged_volume / body_volume, body_density / fluid_density)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(body_volume_=body_volume, body_density_=body_density, fluid_density_=fluid_density)
@validate_output(submerged_volume)
def calculate_submerged_volume(body_volume_: Quantity, body_density_: Quantity,
    fluid_density_: Quantity) -> Quantity:
    if convert_to_si(body_density_) > convert_to_si(fluid_density_):
        raise ValueError("body density must be no greater than the fluid density")

    result_expr = solve(law, submerged_volume)[0]
    result_volume = result_expr.subs({
        body_volume: body_volume_,
        body_density: body_density_,
        fluid_density: fluid_density_,
    })
    return Quantity(result_volume)

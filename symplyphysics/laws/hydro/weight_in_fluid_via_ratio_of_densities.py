"""
Weight in fluid via ratio of densities
======================================

The *Archimedean force* acting on a body immersed in a fluid is equal to the weight of the fluid displaced by the body.
It can be derived that the weight of the body immersed in the fluid is proportional to its weight in vacuum and also depends
on the ratio of the fluid density and body density.

**Links:**

#. `Physics LibreTexts, derivable from here <https://phys.libretexts.org/Bookshelves/University_Physics/Physics_(Boundless)/10%3A_Fluids/10.3%3A_Archimedes_Principle>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    clone_as_symbol,
    symbols,
    Quantity,
    validate_input,
    validate_output,
)

weight_in_fluid = clone_as_symbol(symbols.force,
    display_symbol="W_fl",
    display_latex="W_\\text{fl}")
"""
Weight of the body immersed in the fluid. See :symbols:`force`.
"""

weight_in_vacuum = clone_as_symbol(symbols.force,
    display_symbol="W_vac",
    display_latex="W_\\text{vac}")
"""
Weight of the body in vacuum, i.e. its true weight. See :symbols:`force`.
"""

fluid_density = clone_as_symbol(symbols.density, display_symbol="rho_fl", display_latex="\\rho_\\text{fl}")
"""
:symbols:`density` of the fluid.
"""

body_density = clone_as_symbol(symbols.density, display_symbol="rho_b", display_latex="\\rho_\\text{b}")
"""
:symbols:`density` of the body.
"""

law = Eq(weight_in_fluid, weight_in_vacuum * (1 - (fluid_density / body_density)))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(weight_air_=weight_in_vacuum,
    liquid_density_=fluid_density,
    body_density_=body_density)
@validate_output(weight_in_fluid)
def calculate_weight(weight_air_: Quantity, liquid_density_: Quantity,
    body_density_: Quantity) -> Quantity:
    result_expr = solve(law, weight_in_fluid, dict=True)[0][weight_in_fluid]
    result_weight = result_expr.subs({
        weight_in_vacuum: weight_air_,
        fluid_density: liquid_density_,
        body_density: body_density_,
    })

    return Quantity(result_weight)

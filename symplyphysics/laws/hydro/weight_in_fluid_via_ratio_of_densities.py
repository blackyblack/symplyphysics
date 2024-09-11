"""
Weight in fluid via ratio of densities
======================================

The *Archimedean force* acting on a body immersed in a fluid is equal to the weight of the fluid displaced by the body.
It can be derived that the weight of the body immersed in the fluid is proportional to its weight in vacuum and also depends
on the ratio of the fluid density and body density.
"""

from sympy import (Eq, solve)
from symplyphysics import (clone_symbol, symbols, units, Quantity, Symbol, validate_input,
    validate_output)

weight_in_fluid = clone_symbol(symbols.force,
    display_symbol="W_fl",
    display_latex="W_\\text{fl}")
"""
Weight of the body immersed in the fluid.
"""

weight_in_vacuum = clone_symbol(symbols.force,
    display_symbol="W_vac",
    display_latex="W_\\text{vac}")
"""
Weight of the body in vacuum, i.e. its true weight.
"""

fluid_density = Symbol("fluid_density", units.mass / units.volume)
r"""
Density of the fluid.

Symbol:
    :code:`rho_fl`

Latex:
    :math:`\rho_\text{fl}`
"""

body_density = Symbol("body_density", units.mass / units.volume)
r"""
Density of the body.

Symbol:
    :code:`rho_b`

Latex:
    :math:`\rho_\text{b}`
"""

law = Eq(weight_in_fluid, weight_in_vacuum * (1 - (fluid_density / body_density)))
r"""
:code:`W_fl = W_vac * (1 - rho_fl / rho_b)`

Latex:
    .. math::
        W_\text{fl} = W_\text{vac} \left( 1 - \frac{\rho_\text{fl}}{\rho_\text{b}} \right)
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

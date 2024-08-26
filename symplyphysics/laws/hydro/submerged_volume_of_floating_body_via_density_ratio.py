r"""
Submerged volume of floating body via density ratio
===================================================

If a body is fully or partially submerged in a fluid, an Archimedes force
starts acting on it from the fluid pushing it in the opposite direction of
the gravity force.

The ratio of the body's volume submerged in a fluid to its total volume
depends on the ratio of the densities of the body and the fluid.

**Conditions:**

#. :math:`\rho \le \rho_\text{fl}`, so the body must be floating.
#. The body must be in static equilibrium.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output)
from symplyphysics.core.convert import convert_to_si

submerged_volume = Symbol("submerged_volume", units.volume)
r"""
Volume submerged in the fluid, which is equal to the volume of the displaced fluid.

Symbol:
    :code:`V_fl`

Latex:
    :math:`V_\text{fl}`
"""

body_volume = Symbol("body_volume", units.volume)
"""
Total volume of the body.

Symbol:
    :code:`V`
"""

body_density = Symbol("body_density", units.mass / units.volume)
r"""
Density of the body.

Symbol:
    :code:`rho`

Latex:
    :math:`\rho`
"""

fluid_density = Symbol("fluid_density", units.mass / units.volume)
r"""
Density of the fluid.

Symbol:
    :code:`rho_fl`

Latex:
    :math:`\rho_\text{fl}`
"""

law = Eq(submerged_volume / body_volume, body_density / fluid_density)
r"""
:code:`V_fl / V = rho / rho_fl`

Latex:
    .. math::
        \frac{V_\text{fl}}{V} = \frac{\rho}{\rho_\text{fl}}
"""


@validate_input(
    body_volume_=body_volume,
    body_density_=body_density,
    fluid_density_=fluid_density)
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

r"""
Electric field of uniformly charged plane
=========================================

The electric field strength of a uniformly charged plane is proportional to its
charge density.

**Notation:**

#. :math:`\varepsilon_0` (:code:`epsilon_0`) is vacuum permittivity.

**Conditions:**

#. The plane is thin, i.e. its thickness approaches zero.
"""

from sympy import (Eq, solve)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

electric_field_strength = Symbol("electric_field_strength", units.voltage / units.length)
"""
Value of the electric field.

Symbol:
    :code:`E`
"""

surface_charge_density = Symbol("surface_charge_density", units.charge / units.area)
r"""
Surface charge density of the plane.

Symbol:
    :code:`sigma`

Latex:
    :math:`\sigma`
"""

law = Eq(electric_field_strength, surface_charge_density / (2 * units.vacuum_permittivity))
r"""
:code:`E = sigma / (2 * epsilon_0)`

Latex:
    .. math::
        E = \frac{\sigma}{2 \varepsilon_0}
"""


@validate_input(surface_charge_density_=surface_charge_density)
@validate_output(electric_field_strength)
def calculate_electric_intensity(surface_charge_density_: Quantity) -> Quantity:
    result_expr = solve(law, electric_field_strength, dict=True)[0][electric_field_strength]
    result_expr = result_expr.subs({
        surface_charge_density: surface_charge_density_,
    })
    return Quantity(result_expr)

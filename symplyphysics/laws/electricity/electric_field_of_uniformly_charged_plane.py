r"""
Electric field of uniformly charged plane
=========================================

The electric field strength of a uniformly charged plane is proportional to its
charge density.

**Notation:**

#. :math:`\varepsilon_0` (:code:`epsilon_0`) is vacuum permittivity.

**Conditions:**

#. The plane is thin, i.e. its thickness approaches zero.

**Links:**

#. `Electric field of a uniformly charged plane <https://farside.ph.utexas.edu/teaching/316/lectures/node27.html>`__
"""

from sympy import (Eq, solve)
from symplyphysics import (
    units,
    quantities,
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
)
from symplyphysics.laws.electricity import (
    electric_flux_through_closed_surface_via_total_charge as _gauss_law,
)

electric_field_strength = SymbolNew("E", units.voltage / units.length)
"""
Value of the electric field.
"""

surface_charge_density = SymbolNew("sigma", units.charge / units.area, display_latex="\\sigma")
"""
Surface charge density of the plane.
"""

law = Eq(electric_field_strength, surface_charge_density / (2 * quantities.vacuum_permittivity))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(surface_charge_density_=surface_charge_density)
@validate_output(electric_field_strength)
def calculate_electric_intensity(surface_charge_density_: Quantity) -> Quantity:
    result_expr = solve(law, electric_field_strength, dict=True)[0][electric_field_strength]
    result_expr = result_expr.subs({
        surface_charge_density: surface_charge_density_,
    })
    return Quantity(result_expr)

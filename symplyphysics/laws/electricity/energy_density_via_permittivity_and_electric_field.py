r"""
Energy density via permittivity and electric field
==================================================

Volumetric energy density of the electric field depends on the permittivity of the medium and
on the strength (value) of the electric field at that point.

**Notation:**

#. :math:`\varepsilon_0` (:code:`epsilon_0`) is vacuum permittivity.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output, dimensionless)

energy_density = Symbol("energy_density", units.energy / units.volume)
"""
Volumetric energy density of the electric field, see :doc:`laws.quantities.quantity_is_volumetric_density_times_volume`.

Symbol:
    :code:`w`
"""

relative_permittivity = Symbol("relative_permittivity", dimensionless)
r"""
Relative permittivity of the medium.

Symbol:
    :code:`epsilon`

Latex:
    :math:`\varepsilon`
"""

electric_field = Symbol("electric_field", units.voltage / units.length)
"""
Value of the electric field.

Symbol:
    :code:`E`
"""

law = Eq(energy_density, (units.vacuum_permittivity * relative_permittivity * electric_field**2) / 2)
r"""
:code:`w = 1/2 * epsilon_0 * epsilon * E^2`

Latex:
    .. math::
        w = \frac{1}{2} \varepsilon_0 \varepsilon E^2
"""


@validate_input(relative_permittivity_=relative_permittivity,
    electric_intensity_=electric_field)
@validate_output(energy_density)
def calculate_energy_density(relative_permittivity_: float,
    electric_intensity_: Quantity) -> Quantity:
    result_expr = solve(law, energy_density, dict=True)[0][energy_density]
    result_expr = result_expr.subs({
        relative_permittivity: relative_permittivity_,
        electric_field: electric_intensity_,
    })
    return Quantity(result_expr)
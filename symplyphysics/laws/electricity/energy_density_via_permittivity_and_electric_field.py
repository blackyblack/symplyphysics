r"""
Energy density via permittivity and electric field
==================================================

Volumetric energy density of the electric field depends on the permittivity of the medium and
on the strength of the electric field at that point.

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

absolute_permittivity = Symbol("absolute_permittivity", units.capacitance / units.length)
r"""
Relative permittivity of the medium.

Symbol:
    :code:`epsilon`

Latex:
    :math:`\varepsilon`
"""

electric_field_strength = Symbol("electric_field_strength", units.voltage / units.length)
"""
Strength of the electric field.

Symbol:
    :code:`E`
"""

law = Eq(energy_density, (absolute_permittivity * electric_field_strength**2) / 2)
r"""
:code:`w = 1/2 * epsilon * E^2`

Latex:
    .. math::
        w = \frac{1}{2} \varepsilon E^2
"""


@validate_input(absolute_permittivity_=absolute_permittivity,
    electric_intensity_=electric_field_strength)
@validate_output(energy_density)
def calculate_energy_density(absolute_permittivity_: float,
    electric_intensity_: Quantity) -> Quantity:
    result_expr = solve(law, energy_density, dict=True)[0][energy_density]
    result_expr = result_expr.subs({
        absolute_permittivity: absolute_permittivity_,
        electric_field_strength: electric_intensity_,
    })
    return Quantity(result_expr)

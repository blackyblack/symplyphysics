r"""
Energy density via permittivity and electric field
==================================================

Volumetric energy density of the electric field depends on the permittivity of the medium and
on the strength of the electric field at that point.
"""

from sympy import (Eq, solve)
from symplyphysics import (Quantity, validate_input,
    validate_output, symbols)

energy_density = symbols.energy_density
"""
Volumetric :symbols:`energy_density` of the electric field, see
:doc:`laws.quantities.quantity_is_volumetric_density_times_volume`.
"""

absolute_permittivity = symbols.absolute_permittivity
"""
:symbols:`absolute_permittivity` of the medium.
"""

electric_field_strength = symbols.electric_field_strength
"""
:symbols:`electric_field_strength`.
"""

law = Eq(energy_density, (absolute_permittivity * electric_field_strength**2) / 2)
r"""
:laws:symbol::

:laws:latex::
"""


@validate_input(absolute_permittivity_=absolute_permittivity,
    electric_intensity_=electric_field_strength)
@validate_output(energy_density)
def calculate_energy_density(absolute_permittivity_: Quantity,
    electric_intensity_: Quantity) -> Quantity:
    result_expr = solve(law, energy_density, dict=True)[0][energy_density]
    result_expr = result_expr.subs({
        absolute_permittivity: absolute_permittivity_,
        electric_field_strength: electric_intensity_,
    })
    return Quantity(result_expr)

"""
Refractive index via permittivity and permeability
==================================================

Refractive index can be calculated from the relative permittivity and permeability
of the medium.
"""

from sympy import Eq, solve, sqrt
from symplyphysics import validate_input, validate_output, convert_to_float, symbols

refractive_index = symbols.relative_refractive_index
"""
:symbols:`relative_refractive_index` of the medium.
"""

relative_permittivity = symbols.relative_permittivity
"""
:symbols:`relative_permittivity` of the medium.
"""

relative_permeability = symbols.relative_permeability
"""
:symbols:`relative_permeability` of the medium.
"""

law = Eq(refractive_index, sqrt(relative_permittivity * relative_permeability))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    relative_dielectric_permeability_=relative_permittivity,
    relative_magnetic_permeability_=relative_permeability,
)
@validate_output(refractive_index)
def calculate_refraction_factor(relative_dielectric_permeability_: float,
    relative_magnetic_permeability_: float) -> float:
    result_expr = solve(law, refractive_index, dict=True)[0][refractive_index]
    factor_applied = result_expr.subs({
        relative_permittivity: relative_dielectric_permeability_,
        relative_permeability: relative_magnetic_permeability_
    })
    return convert_to_float(factor_applied)

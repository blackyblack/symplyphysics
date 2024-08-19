"""
Refractive index via permittivity and permeability
==================================================

Refractive index can be calculated from the relative permittivity and permeability
of the medium.
"""

from sympy import (Eq, solve, sqrt)
from symplyphysics import dimensionless, Symbol, validate_input, validate_output, convert_to_float

refractive_index = Symbol("refractive_index", dimensionless)
"""
Refractive index of the medium.

Symbol:
    :code:`n`
"""

relative_permittivity = Symbol("relative_permittivity", dimensionless)
r"""
Relative permittivity of the medium.

Symbol:
    :code:`epsilon_r`

Latex:
    :math:`\varepsilon_\text{r}`
"""

relative_permeability = Symbol("relative_permeability", dimensionless)
r"""
Relative permeability of the medium.

Symbol:
    :code:`mu_r`

Latex:
    :math:`\mu_\text{r}`
"""

law = Eq(refractive_index, sqrt(relative_permittivity * relative_permeability))
r"""
:code:`n = sqrt(epsilon_r * mu_r)`

Latex:
    .. math::
        n = \sqrt{\varepsilon_\text{r} \mu_\text{r}}
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

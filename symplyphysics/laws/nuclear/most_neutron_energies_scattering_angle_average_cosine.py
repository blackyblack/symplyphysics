"""
Average cosine of scattering angle from mass number
===================================================

Averaged cosine of the angle of neutron scattering depends on the mass number of the
target nucleus.

**Links:**

#. `NuclearPower <https://www.nuclear-power.com/nuclear-power/reactor-physics/neutron-diffusion-theory/diffusion-coefficient/>`__.
"""

from sympy import (Eq, solve)
from symplyphysics import symbols, Symbol, dimensionless

mass_number = symbols.mass_number
"""
:symbols:`mass_number` of the target nucleus.
"""

average_scattering_angle_cosine = Symbol("mu", dimensionless, display_latex="\\mu")
"""
Average of the cosine of the angle at which neutrons are scattered in the medium in the
lab system.
"""

law = Eq(average_scattering_angle_cosine, 2 / (3 * mass_number))
"""
:laws:symbol::

:laws:latex::
"""


def calculate_average_scattering_angle_cosine(target_nucleus_mass_number_: int) -> float:
    result_angle_cosine_expr = solve(law, average_scattering_angle_cosine,
        dict=True)[0][average_scattering_angle_cosine]
    result_expr = result_angle_cosine_expr.subs(mass_number,
        target_nucleus_mass_number_)
    return float(result_expr.evalf())

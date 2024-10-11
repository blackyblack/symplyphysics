"""
Wave speed from medium permittivity and permeability
====================================================

Speed of electromagnetic wave propagation depends on relative permittivity and relative
permeability of the medium.

**Notation:**

#. :quantity_notation:`speed_of_light`.
"""

from sympy import (Eq, solve, sqrt)
from symplyphysics import (units, Quantity, Symbol, dimensionless, validate_input, validate_output,
    quantities)

wave_speed = Symbol("wave_speed", units.velocity)
"""
Speed of wave propagation in medium.

Symbol:
    :code:`v`
"""

relative_permittivity = Symbol("relative_permittivity", dimensionless)
r"""
Permittivity of medium relative to that of vacuum.

Symbol:
    :code:`epsilon`

Latex:
    :math:`\varepsilon`
"""

relative_permeability = Symbol("relative_permeability", dimensionless)
r"""
Permeability of medium relative to that of vacuum.

Symbol:
    :code:`mu`

Latex:
    :math:`\mu`
"""

law = Eq(wave_speed,
    quantities.speed_of_light / sqrt(relative_permittivity * relative_permeability))
r"""
:code:`v = c / sqrt(epsilon * mu)`

Latex:
    .. math::
        v = \frac{c}{\sqrt {\varepsilon \mu}}
"""


@validate_input(permittivity_=relative_permittivity, permeability_=relative_permeability)
@validate_output(wave_speed)
def calculate_wavespeed(permittivity_: float, permeability_: float) -> Quantity:
    result_expr = solve(law, wave_speed, dict=True)[0][wave_speed]
    wavespeed_applied = result_expr.subs({
        relative_permittivity: permittivity_,
        relative_permeability: permeability_
    })
    return Quantity(wavespeed_applied)

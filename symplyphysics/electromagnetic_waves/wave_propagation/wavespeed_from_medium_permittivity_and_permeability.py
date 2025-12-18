"""
Wave speed from medium permittivity and permeability
====================================================

Speed of electromagnetic wave propagation depends on relative permittivity and relative
permeability of the medium.

**Notation:**

#. :quantity_notation:`speed_of_light`.

**Links:**

#. `Wikipedia, derivable from here <https://en.wikipedia.org/wiki/Phase_velocity#Refractive_index>`__.
"""

from sympy import (Eq, solve, sqrt)
from symplyphysics import (Quantity, validate_input, validate_output, quantities, symbols)

wave_speed = symbols.phase_speed
"""
Speed of wave propagation in medium. See :symbols:`phase_speed`.
"""

relative_permittivity = symbols.relative_permittivity
"""
:symbols:`relative_permittivity` of the medium.
"""

relative_permeability = symbols.relative_permeability
"""
:symbols:`relative_permeability` of the medium.
"""

law = Eq(wave_speed,
    quantities.speed_of_light / sqrt(relative_permittivity * relative_permeability))
"""
:laws:symbol::

:laws:latex::
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

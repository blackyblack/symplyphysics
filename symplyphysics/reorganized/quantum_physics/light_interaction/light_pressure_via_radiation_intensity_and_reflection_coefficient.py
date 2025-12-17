"""
Light pressure
==============

Light pressure is the pressure exerted by light radiation incident on the surface of a body.
Under normal light incidence on the surface, the light pressure depends only on the reflection
coefficient and the intensity of the radiation.

**Notation:**

#. :quantity_notation:`speed_of_light`.

**Links:**

#. `BYJU's, "Radiation pressure formula" <https://byjus.com/physics/radiation-pressure/>`__.
"""

from sympy import (Eq, solve)
from symplyphysics import (Quantity, validate_input, validate_output, symbols, quantities)

pressure = symbols.pressure
"""
Light :symbols:`pressure`.
"""

intensity = symbols.intensity
"""
:symbols:`intensity` of incident radiation.
"""

reflection_coefficient = symbols.reflectance
"""
:symbols:`reflectance` of the surface.
"""

law = Eq(pressure, intensity * (1 + reflection_coefficient) / quantities.speed_of_light)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(intensity_=intensity, reflection_coefficient_=reflection_coefficient)
@validate_output(pressure)
def calculate_pressure(intensity_: Quantity, reflection_coefficient_: float) -> Quantity:
    result_expr = solve(law, pressure, dict=True)[0][pressure]
    result_expr = result_expr.subs({
        intensity: intensity_,
        reflection_coefficient: reflection_coefficient_,
    })
    return Quantity(result_expr)

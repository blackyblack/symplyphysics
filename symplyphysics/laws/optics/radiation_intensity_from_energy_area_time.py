"""
Radiation intensity from energy, area, and time
===============================================

Intensity is a scalar physical quantity that quantitatively characterizes the power carried by
a wave in the direction of propagation per unit area per unit time.

**Links:**

#. `Wikipedia, derivable from here <https://en.wikipedia.org/wiki/Intensity_(physics)#Mathematical_description>`__.
"""

from sympy import Eq, solve
from symplyphysics import Quantity, validate_input, validate_output, symbols

intensity = symbols.intensity
"""
(Average) :symbols:`intensity` of incident radiation.
"""

energy = symbols.energy
"""
:symbols:`energy` of incident radiation.
"""

area = symbols.area
"""
:symbols:`area`.
"""

time = symbols.time
"""
:symbols:`time`.
"""

law = Eq(intensity, energy / (area * time))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    energy_=energy,
    area_=area,
    time_=time,
)
@validate_output(intensity)
def calculate_intensity(energy_: Quantity, area_: Quantity, time_: Quantity) -> Quantity:
    result_expr = solve(law, intensity, dict=True)[0][intensity]
    intensity_applied = result_expr.subs({energy: energy_, area: area_, time: time_})
    return Quantity(intensity_applied)

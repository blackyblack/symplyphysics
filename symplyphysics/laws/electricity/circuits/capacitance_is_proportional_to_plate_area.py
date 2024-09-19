r"""
Capacitance is proportional to plate area
=========================================

If a capacitor is composed of two equal parallel plates, its capacitance is proportional to
the area of the plates and the inverse of the distance between the plates. It also depends
on the permittivity of the medium that fills the space between the plates.

**Conditions:**

#. The plates must be large enough for the electric field to be constant between the plates.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

capacitance = symbols.capacitance
"""
:symbols:`capacitance` of the capacitor.
"""

absolute_permittivity = symbols.absolute_permittivity
"""
:symbols:`absolute_permittivity` of the medium between the plates.
"""

area = symbols.area
"""
:symbols:`area` of the plates.
"""

distance = symbols.distance
"""
:symbols:`distance` between the plates.
"""

law = Eq(capacitance, absolute_permittivity * area / distance)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(plate_area_=area, distance_between_plates_=distance)
@validate_output(capacitance)
def calculate_capacitance(dielectric_permeability_: float, plate_area_: Quantity,
    distance_between_plates_: Quantity) -> Quantity:
    result_capacitance_expr = solve(law, capacitance, dict=True)[0][capacitance]
    result_expr = result_capacitance_expr.subs({
        absolute_permittivity: dielectric_permeability_,
        area: plate_area_,
        distance: distance_between_plates_
    })
    return Quantity(result_expr)

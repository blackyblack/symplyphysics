r"""
Capacitance is proportional to plate area
=========================================

If a capacitor is composed of two equal parallel plates, its capacitance is proportional to
the area of the plates and the inverse of the distance between the plates.

**Conditions:**

#. The plates must be large enough for the electric field to be constant between the plates.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output)

capacitance = Symbol("capacitance", units.capacitance)
"""
Capacitance of the capacitor.

Symbol:
    :code:`C`
"""

absolute_permittivity = Symbol("absolute_permittivity", units.capacitance / units.length)
r"""
Absolute permittivity of the medium between the plates.

Symbol:
    :code:`epsilon`

Latex:
    :math:`\varepsilon`
"""

area = Symbol("area", units.area)
"""
Area of the plates.

Symbol:
    :code:`A`
"""

distance = Symbol("distance", units.length)
"""
Distance between the plates.

Symbol:
    :code:`d`
"""

law = Eq(capacitance,
    absolute_permittivity * area / distance)
r"""
:code:`C = epsilon * A / d`

Latex:
    .. math::
        C = \frac{\varepsilon A}{d}
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

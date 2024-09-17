r"""
Electric field outside charged sphere
=====================================

The electric field outside of a charged sphere behaves as is the sphere were a point charge,
i.e. its magnitude is inversely proportional to the square of the distance to the center
of the sphere. However, due to the Gauss's law, on the inside the electric field is exactly zero.

**Notation:**

#. :math:`\varepsilon_0` (:code:`epsilon_0`) is vacuum permittivity.

**Conditions:**

#. The sphere is thin, i.e. its thickness approaches zero.
"""

from sympy import (Eq, solve, pi)
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    quantities,
)

electric_field_strength = symbols.electric_field_strength
"""
Value of the electric field.
"""

charge = symbols.charge
"""
Total charge of the sphere.
"""

distance = symbols.distance_to_origin
"""
Distance to the center of the sphere.
"""

law = Eq(electric_field_strength, charge / (4 * pi * quantities.vacuum_permittivity * distance**2))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(charge_=charge, distance_=distance)
@validate_output(electric_field_strength)
def calculate_electric_intensity(charge_: Quantity, distance_: Quantity) -> Quantity:
    result_expr = solve(law, electric_field_strength, dict=True)[0][electric_field_strength]
    result_expr = result_expr.subs({
        charge: charge_,
        distance: distance_,
    })
    return Quantity(result_expr)

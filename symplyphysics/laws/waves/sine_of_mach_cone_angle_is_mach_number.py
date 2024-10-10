r"""
Sine of Mach cone angle via Mach number
=======================================

If the speed of a source relative to the medium exceeds the speed of sound in the medium,
the Doppler equation no longer applies and this results in shock waves. The wavefronts
of the waves originating from the source form a cone, namely a *Mach cone*. The half-angle
of the cone is called the *Mach cone angle*, which is related to the Mach number of the source.
See the `illustration <https://www.grc.nasa.gov/www/k-12/airplane/machang.html>`_ of the
phenomenon.

**Conditions:**

#. :math:`M \ge 1`, i.e. the source speed exceeds the speed of sound in the medium.
"""

from sympy import Eq, sin, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

mach_cone_angle = symbols.angle
"""
Mach cone :symbols:`angle`, which is the angle between the Mach wave wavefront (the Mach cone) and
the vector pointing opposite to the velocity vector of the source.
"""

mach_number = symbols.mach_number
"""
:symbols:`mach_number` of the moving source.
"""

law = Eq(sin(mach_cone_angle), 1 / mach_number)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(mach_number_=mach_number)
@validate_output(mach_cone_angle)
def calculate_mach_cone_angle(mach_number_: float) -> Quantity:
    if mach_number_ < 1:
        raise ValueError("The Mach number must be greater or equal to 1")

    result_expr = solve(law, mach_cone_angle)[1]
    result = result_expr.subs(mach_number, mach_number_)
    return Quantity(result)

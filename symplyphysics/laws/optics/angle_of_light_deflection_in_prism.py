"""
Angle of light deflection in prism
==================================

A prism, an optical prism, is a body made of a homogeneous material transparent to optical radiation,
bounded by flat reflecting and refractive surfaces located at strictly defined angles to each other.
With a small angle of incidence of the incoming beam, the angle of deflection of the beam depends only
on the angle between the faces of the prism and the refractive index of the prism.
There is a ray falling on the prism and a ray coming out of the prism. Let's continue these rays inside
the prism. Then the angle of intersection of these rays, looking towards the output beam, will be called
the angle of deviation.

**Conditions:**

#. The deviation angle and the incident angle are small.

**Links:**

#. `Wikipedia <https://ru.wikipedia.org/wiki/Призма_(оптика)#:~:text=Призма%2C%20оптическая%20призма%20—%20тело%20из,определёнными%20углами%20друг%20к%20другу>`__.

..
    TODO find link in English
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

deviation_angle = clone_as_symbol(symbols.angle, display_symbol="b", display_latex="b")
"""
:symbols:`angle` deviation in the prism.
"""

face_angle = clone_as_symbol(symbols.angle, display_symbol="a", display_latex="a")
"""
:symbols:`angle` between faces of the prism.
"""

relative_refractive_index = symbols.relative_refractive_index
"""
:symbols:`relative_refractive_index` of the prism.
"""

law = Eq(deviation_angle, (face_angle * (relative_refractive_index - 1)))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(angle_faces_=face_angle, refractive_index_=relative_refractive_index)
@validate_output(deviation_angle)
def calculate_angle_deviation(angle_faces_: float | Quantity, refractive_index_: float) -> Quantity:
    result_expr = solve(law, deviation_angle, dict=True)[0][deviation_angle]
    result_expr = result_expr.subs({
        face_angle: angle_faces_,
        relative_refractive_index: refractive_index_,
    })
    return Quantity(result_expr)

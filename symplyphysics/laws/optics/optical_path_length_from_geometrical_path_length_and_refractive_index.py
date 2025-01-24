"""
Optical path length from geometrical path length and refractive index
=====================================================================

*Geometrical path length* is the Euclidean distance integrated along a light ray
between any two points.

**Conditions:**

#. The environment must be homogeneous.

**Links:**

#. `Wikipedia, first formula <https://en.wikipedia.org/wiki/Optical_path_length#>`__.
"""

from sympy import (Eq, solve)
from symplyphysics import (Quantity, validate_input, validate_output, symbols)

optical_path = symbols.optical_distance
"""
Optical path length of the light ray. See :symbols:`optical_distance`.
"""

geometrical_path = symbols.distance
"""
Geometrical path length of the light ray. See :symbols:`distance`.
"""

refractive_index = symbols.relative_refractive_index
"""
:symbols:`relative_refractive_index` of the medium.
"""

law = Eq(optical_path, refractive_index * geometrical_path)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(geometric_path_=geometrical_path, refractive_index_=refractive_index)
@validate_output(optical_path)
def calculate_optical_path(geometric_path_: Quantity, refractive_index_: float) -> Quantity:
    result_expr = solve(law, optical_path, dict=True)[0][optical_path]
    result_optical_path = result_expr.subs({
        geometrical_path: geometric_path_,
        refractive_index: refractive_index_,
    })
    return Quantity(result_optical_path)

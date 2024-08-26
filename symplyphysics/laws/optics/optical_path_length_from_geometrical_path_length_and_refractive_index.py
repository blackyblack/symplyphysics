"""
Optical path length from geometrical path length and refractive index
=====================================================================

*Optical path length*, also known as *optical length* or *optical distance*, is the
length that light needs to travel through a vacuum to create the same phase difference
as it would have when traveling through a given medium.

*Geometrical path length* is the Euclidean distance integrated along a light ray
between any two points.

**Conditions:**

#. The environment must be homogeneous.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output, dimensionless)

optical_path = Symbol("optical_path", units.length)
"""
Optical path length of the light ray.

Symbol:
    :code:`L`
"""

geometrical_path = Symbol("geometrical_path", units.length)
"""
Geometrical path length of the light ray.

Symbol:
    :code:`l`
"""

refractive_index = Symbol("refractive_index", dimensionless)
"""
Refractive index of the medium.

Symbol:
    :code:`n`
"""

law = Eq(optical_path, refractive_index * geometrical_path)
r"""
:code:`L = n * l`

Latex:
    .. math::
        L = n l
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

"""
Divergence of magnetic flux density is zero
===========================================

The divergence of magnetic induction is equal to zero. Another form of this law states that
magnetic charges do not exist.

**Links:**

#. `Wikipedia, second equation <https://en.wikipedia.org/wiki/Maxwell%27s_equations#>`__.
"""

from sympy import Eq
from symplyphysics import symbols

from symplyphysics.core.experimental.vectors import (clone_as_vector_function,
    clone_as_vector_symbol)
from symplyphysics.core.experimental.operators import VectorDivergence

position_vector = clone_as_vector_symbol(symbols.distance_to_origin)
"""
Position vector of the point in space. See :symbols:`distance_to_origin`.
"""

magnetic_flux_density = clone_as_vector_function(
    symbols.magnetic_flux_density,
    (position_vector,),
)
"""
Vector field of the :symbols:`magnetic_flux_density` as a function of the
:attr:`~position_vector`.
"""

law = Eq(
    VectorDivergence(magnetic_flux_density(position_vector), evaluate=False),
    0,
)
"""
..
    NOTE: code printers have not been implemented yet for `VectorDivergence`

:code:`div(B(r)) = 0`

Latex:
    .. math::
        \\text{div} \\, {\\vec B} \\left( \\vec r \\right) = 0
"""

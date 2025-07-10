"""
Mechanical work is line integral of force
=========================================

Mechanical work is a scalar quantity denoting the energy transferred to or from an object via the
application of force along a displacement.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Work_(physics)#Mathematical_calculation>`__.
"""

from sympy import Eq, Expr
from symplyphysics import Quantity, validate_output, symbols, Symbol
from symplyphysics.core.symbols.symbols import BasicSymbol
from symplyphysics.core.dimensions import any_dimension

from symplyphysics.core.experimental.vectors import (clone_as_vector_function,
    clone_as_vector_symbol, VectorDot)
from symplyphysics.core.experimental.coordinate_systems.curve import Curve
from symplyphysics.core.experimental.integrals.line_integral import (LineIntegral,
    INFINITESIMAL_DISPLACEMENT as dr)

work = symbols.work
"""
Mechanical :symbols:`work` done by the :attr:`~force` to displace the body.
"""

position_vector = clone_as_vector_symbol(symbols.distance_to_origin)
"""
Position vector of the body. See :symbols:`distance_to_origin`.
"""

force = clone_as_vector_function(symbols.force, (position_vector,))
"""
Vector of the exerted force as a function of :attr:`~position_vector`.
"""

curve = BasicSymbol("C")
"""
Curve given by the body's trajectory during its displacement.
"""

initial_parameter = Symbol("t_0", any_dimension)
"""
Initial value of the curve parameter.
"""

final_parameter = Symbol("t_0", any_dimension)
"""
Final value of the curve parameter.
"""

law = Eq(
    work,
    LineIntegral(VectorDot(force(position_vector), dr), curve,
    (initial_parameter, final_parameter)),
)
"""
:laws:symbol::

:laws:latex::
"""

# Derivable from `./mechanical_work_from_force_and_move` by summing over a piecewise approximation
# of the curve and taking the limit when the length of the pieces approaches zero.


@validate_output(work)
def calculate_work(
    force_: Expr,
    curve_: Curve,
    initial_parameter_: Quantity,
    final_parameter_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        force(position_vector): force_,
        curve: curve_,
        initial_parameter: initial_parameter_,
        final_parameter: final_parameter_,
    }).doit().doit()

    return Quantity(result)

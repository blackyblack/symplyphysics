"""
Circulation is integral along curve
===================================

The circulation of a vector field along a closed curve is defined via a curvilinear integral of
the field along the curve.

**Conditions:**

#. The parametrization of the curve is a vector function of one parameter that is continuously
   differentiable with respect to it.
#. The curve is closed.
"""

from sympy import Eq
from symplyphysics import Quantity
from symplyphysics.core.dimensions import any_dimension, assert_equivalent_dimension
from symplyphysics.core.symbols.symbols import BasicSymbol, Symbol
from symplyphysics.core.coordinate_systems import CoordinateVector
from symplyphysics.core.coordinate_systems.curve import Curve
from symplyphysics.core.vectors import VectorSymbol, VectorDot
from symplyphysics.core.integrals.line_integral import LineIntegral, INFINITESIMAL_DISPLACEMENT

circulation = Symbol("G", any_dimension)
"""
Circulation of the vector :attr:`~field`.
"""

field = VectorSymbol("F", dimension=any_dimension)
"""
Any vector field, i.e. a vector-valued function that depends on the position vector.
"""

curve = BasicSymbol("C")
"""
Curve along which the :attr:`~circulation` is calculated.
"""

initial_parameter = Symbol("u_1", any_dimension)
"""
Initial value of the curve parameter.
"""

final_parameter = Symbol("u_2", any_dimension)
"""
Final value of the curve parameter.
"""

law = Eq(
    circulation,
    LineIntegral(VectorDot(field, INFINITESIMAL_DISPLACEMENT), curve,
    (initial_parameter, final_parameter)),
)
"""
:laws:symbol::

:laws:latex::
"""


def calculate_circulation(
    field_: CoordinateVector,
    curve_: Curve,
    bounds_: tuple[Quantity, Quantity],
) -> Quantity:
    u1_, u2_ = bounds_
    assert_equivalent_dimension(u1_, "initial_parameter_", "calculate_circulation", u2_)

    result = law.rhs.subs({
        field: field_,
        curve: curve_,
        initial_parameter: u1_,
        final_parameter: u2_,
    }).doit().doit()
    return Quantity(result)

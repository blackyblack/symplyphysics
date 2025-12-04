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
from symplyphysics import Quantity, validate_input
from symplyphysics.core.dimensions import any_dimension, dimensionless
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

initial_parameter = Symbol("u_1", dimensionless)
"""
Initial value of the curve parameter.
"""

final_parameter = Symbol("u_2", dimensionless)
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


@validate_input(
    initial_parameter_=initial_parameter,
    final_parameter_=final_parameter,
)
def calculate_circulation(
    field_: CoordinateVector,
    curve_: Curve,
    initial_parameter_: Quantity | float,
    final_parameter_: Quantity | float,
) -> Quantity:
    result = law.rhs.subs({
        field: field_,
    }).subs({
        curve: curve_,
        initial_parameter: initial_parameter_,
        final_parameter: final_parameter_,
    }).doit().doit()
    return Quantity(result)

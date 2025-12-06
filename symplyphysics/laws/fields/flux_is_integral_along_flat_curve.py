"""
Flux is integral along flat curve
=================================

The flux of a vector field exiting a boundary flat curve is defined as the line integral of the
component of the field normal to the curve along the curve.

**Conditions:**

#. The normal to the curve is outward (see right-hand rule).
#. The vector field is continuously differentiable everywhere within the region of the surface.
#. The curve is closed and flat.
"""

from sympy import Eq
from symplyphysics import Quantity, validate_input
from symplyphysics.core.dimensions import any_dimension, dimensionless
from symplyphysics.core.symbols.symbols import BasicSymbol, Symbol
from symplyphysics.core.coordinate_systems import CoordinateVector
from symplyphysics.core.coordinate_systems.curve import Curve
from symplyphysics.core.vectors import VectorSymbol, VectorDot
from symplyphysics.core.integrals.line_integral import LineIntegral, INFINITESIMAL_ARC_LENGTH

flux = Symbol("H", any_dimension)
"""
Flux of the vector :attr:`~field`.
"""

field = VectorSymbol("F", any_dimension)
"""
Vector field, i.e. a vector-valued function of the position vector.
"""

unit_normal = VectorSymbol("n", any_dimension)
"""
Unit normal vector to the curve.
"""

curve = BasicSymbol("C")
"""
Curve which is the boundary of the surface along which :attr:`~flux` is calculated.

Symbol:
    :code:`C`

Latex:
    :math:`C`
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
    flux,
    LineIntegral(
    VectorDot(field, unit_normal) * INFINITESIMAL_ARC_LENGTH, curve,
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
def calculate_flux(
    field_: CoordinateVector,
    curve_normal_: CoordinateVector,
    curve_: Curve,
    initial_parameter_: Quantity | float,
    final_parameter_: Quantity | float,
) -> Quantity:
    result = law.rhs.subs({
        field: field_,
        unit_normal: curve_normal_,
    }).subs({
        curve: curve_,
        initial_parameter: initial_parameter_,
        final_parameter: final_parameter_,
    }).simplify().doit()

    return Quantity(result)

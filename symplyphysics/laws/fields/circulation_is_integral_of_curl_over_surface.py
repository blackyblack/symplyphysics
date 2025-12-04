"""
Circulation is integral of curl over surface
============================================

The circulation of a vector field, defined as :ref:`such <Circulation is integral along curve>`,
can be re-written as a double integral of the flux of that field over any surface whose boundary
is the curve along which the circulation is calculated.

**Notes:**

#. This law is also known as the `Stokes' theorem <https://en.wikipedia.org/wiki/Stokes%27_theorem>`__.

**Conditions:**

#. The surface is smooth.
#. The orientation of the surface and the orientation of the curve are related by the right-hand
   rule.
"""

from sympy import Eq
from symplyphysics import Quantity
from symplyphysics.core.dimensions import any_dimension, assert_equivalent_dimension
from symplyphysics.core.symbols.symbols import BasicSymbol, Symbol
from symplyphysics.core.coordinate_systems import CoordinateVector
from symplyphysics.core.coordinate_systems.surface import Surface
from symplyphysics.core.vectors import VectorSymbol, VectorDot
from symplyphysics.core.operators import VectorCurl
from symplyphysics.core.integrals.surface_integral import SurfaceIntegral, INFINITESIMAL_VECTOR_AREA

circulation = Symbol("G", any_dimension)
"""
Circulation of the vector :attr:`~field`.
"""

field = VectorSymbol("F", any_dimension)
"""
Any vector field, i.e. a vector-valued function that depends on the position vector.
"""

surface = BasicSymbol("S")
"""
Surface along which the integral is evaluated.

Symbol:

    :code:`S`

Latex:

    :math:`S`
"""

initial_first_parameter = Symbol("u_1", any_dimension)
"""
Initial value of the first surface parameter.
"""

final_first_parameter = Symbol("u_2", any_dimension)
"""
Final value of the first surface parameter.
"""

initial_second_parameter = Symbol("v_1", any_dimension)
"""
Initial value of the second surface parameter.
"""

final_second_parameter = Symbol("v_2", any_dimension)
"""
Final value of the second surface parameter.
"""

law = Eq(
    circulation,
    SurfaceIntegral(
    VectorDot(VectorCurl(field, evaluate=False), INFINITESIMAL_VECTOR_AREA, evaluate=False),
    surface,
    ((initial_first_parameter, final_first_parameter),
    (initial_second_parameter, final_second_parameter)),
    ),
)
r"""
:laws:symbol::

..
    FIXME Unimplemented latex printer for the surface integral

Latex:

    .. math::

        G = \iint \limits_S \left( \nabla \times \vec F \right) \cdot d \vec S
"""


def calculate_circulation(
    field_: CoordinateVector,
    surface_: Surface,
    bounds_: tuple[tuple[Quantity, Quantity], tuple[Quantity, Quantity]],
) -> Quantity:
    (u1_, u2_), (v1_, v2_) = bounds_
    assert_equivalent_dimension(u1_, "first_initial_parameter_", "calculate_circulation", u2_)
    assert_equivalent_dimension(v1_, "second_initial_parameter_", "calculate_circulation", v2_)
    result = law.rhs.subs({
        field: field_,
    }).subs({
        surface: surface_,
        initial_first_parameter: u1_,
        final_first_parameter: u2_,
        initial_second_parameter: v1_,
        final_second_parameter: v2_,
    })

    result = result.doit().doit()
    return Quantity(result)

"""
Flux is integral of divergence over surface
===========================================

The flux of a vector field exiting a closed flat curve can be re-written as a surface integral of
the divergence of the given field over the surface whose boundary is the given curve.

**Notes:**

#. This is also known as the `Green's theorem <https://en.wikipedia.org/wiki/Green%27s_theorem>`__.

**Conditions:**

#. The vector field is continuously differentiable everywhere within the region of the surface.
#. The curve is closed and flat.
#. The curve and the surface are oriented according to the right-hand rule.
"""

from sympy import Eq
from symplyphysics import Quantity, validate_input
from symplyphysics.core.dimensions import any_dimension, dimensionless
from symplyphysics.core.symbols.symbols import BasicSymbol, Symbol
from symplyphysics.core.coordinate_systems import CoordinateVector, CoordinateScalar
from symplyphysics.core.coordinate_systems.surface import Surface
from symplyphysics.core.vectors import VectorSymbol, VectorNorm
from symplyphysics.core.operators import VectorDivergence
from symplyphysics.core.integrals.surface_integral import SurfaceIntegral, INFINITESIMAL_VECTOR_AREA

flux = Symbol("H", any_dimension)
"""
Flux of the vector :attr:`~field`.
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

initial_first_parameter = Symbol("u_1", dimensionless)
"""
Initial value of the first surface parameter.
"""

final_first_parameter = Symbol("u_2", dimensionless)
"""
Final value of the first surface parameter.
"""

initial_second_parameter = Symbol("v_1", dimensionless)
"""
Initial value of the second surface parameter.
"""

final_second_parameter = Symbol("v_2", dimensionless)
"""
Final value of the second surface parameter.
"""

law = Eq(
    flux,
    SurfaceIntegral(
    VectorDivergence(field, evaluate=False) * VectorNorm(INFINITESIMAL_VECTOR_AREA), surface,
    ((initial_first_parameter, final_first_parameter),
    (initial_second_parameter, final_second_parameter))),
)
r"""
:laws:symbol::

..
    FIXME Implement the surface integral in the latex printer

Latex:

    .. math::

        H = \iint \limits_S \text{div} \vec F \, dS
"""


@validate_input(
    initial_first_parameter_=initial_first_parameter,
    final_first_parameter_=final_first_parameter,
    initial_second_parameter_=initial_second_parameter,
    final_second_parameter_=final_second_parameter,
)
def calculate_flux(
    field_: CoordinateVector,
    surface_: Surface,
    initial_first_parameter_: Quantity | float,
    final_first_parameter_: Quantity | float,
    initial_second_parameter_: Quantity | float,
    final_second_parameter_: Quantity | float,
) -> Quantity:
    result = law.rhs.subs(field, field_).subs({
        surface: surface_,
        initial_first_parameter: initial_first_parameter_,
        final_first_parameter: final_first_parameter_,
        initial_second_parameter: initial_second_parameter_,
        final_second_parameter: final_second_parameter_,
    })

    # Result of applying `VectorDivergence`, we have to get rid of it manually. Otherwise, the
    # integral evaluation below fails.
    for expr in result.atoms(CoordinateScalar):
        result = result.subs(expr, expr.scalar)

    result = result.doit()

    return Quantity(result)

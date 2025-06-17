"""
Curl of electric field is negative magnetic flux density derivative
===================================================================

**Faraday's law of induction** states that a change in magnetic flux density generates a
rotational electric field. This law is valid for any magnetic field that changes over time.

**Links:**

#. `Wikipedia, fourth line in table <https://en.wikipedia.org/wiki/Maxwell%27s_equations#Macroscopic_formulation>`__.

#. `Physics LibreTexts, formula 15.7.1 <https://phys.libretexts.org/Bookshelves/Electricity_and_Magnetism/Electricity_and_Magnetism_(Tatum)/15%3A_Maxwell's_Equations/15.07%3A_Maxwell's_Fourth_Equation>`__.

..
    TODO: rename file
"""

from sympy import Eq, Expr
from symplyphysics import Quantity, validate_input, validate_output, symbols, units
from symplyphysics.core.dimensions.dimensions import assert_quantity_point

from symplyphysics.core.experimental.vectors import (clone_as_vector_function,
    clone_as_vector_symbol, VectorDerivative)
from symplyphysics.core.experimental.operators import VectorCurl
from symplyphysics.core.experimental.coordinate_systems import (QuantityCoordinateVector,
    AppliedPoint, CoordinateVector)
from symplyphysics.core.experimental.solvers import solve_for_vector

time = symbols.time
"""
:symbols:`time`.
"""

position_vector = clone_as_vector_symbol(symbols.distance_to_origin)
"""
Position vector of a point in space. See :symbols:`distance_to_origin`.
"""

electric_field = clone_as_vector_function(
    symbols.electric_field_strength,
    (position_vector, time),
)
"""
Vector of the electric field as a function of :attr:`~position_vector` and :attr:`~time`. See
:symbols:`electric_field_strength`.
"""

magnetic_flux_density = clone_as_vector_function(symbols.magnetic_flux_density,
    (position_vector, time))
"""
Vector of the :symbols:`magnetic_flux_density` field as a function of of :attr:`~position_vector`
and :attr:`~time`.
"""

law = Eq(
    VectorCurl(electric_field(position_vector, time), evaluate=False),
    -1 * VectorDerivative(magnetic_flux_density(position_vector, time), time),
)
"""
..
    NOTE: code printers have not been implemented for `VectorCurl` yet

:code:`curl(E(r, t)) = -1 * Derivative(B(r, t), t)`

Latex:
    .. math::
        \\text{curl} \\, \\vec E \\! \\left( \\vec r, t \\right)
            = - \\frac{\\partial}{\\partial t} \\vec B \\! \\left( \\vec r, t \\right)
"""


@validate_input(time_=time)
@validate_output(units.magnetic_density / units.time)
def calculate_magnetic_induction_derivative_at_point(
    electric_intensity_: Expr,
    point_: AppliedPoint,
    time_: Quantity,
) -> QuantityCoordinateVector:
    assert_quantity_point(point_, "calculate_magnetic_induction_derivative_at_point")

    expr = solve_for_vector(
        law,
        VectorDerivative(magnetic_flux_density(position_vector, time), time),
    )

    result = expr.subs(
        electric_field(position_vector, time),
        electric_intensity_,
    ).doit()
    result = CoordinateVector.from_expr(result)

    result = result.subs(point_.coordinates).subs(time, time_)

    return QuantityCoordinateVector(result.components, result.system, result.point)

"""
Curl of magnetic field is free current density and electric displacement derivative
===================================================================================

The magnetic field circulation theorem states that an electric current and a change in electric
displacement generate a rotational magnetic field. Also known as the **Ampère's circuital law**.

**Notes:**

#. The :math:`\\text{curl}` operator is only defined for a 3D space.

**Links:**

#. `Wikipedia, second line in table <https://en.wikipedia.org/wiki/Maxwell%27s_equations#Macroscopic_formulation>`__.
#. `Physics LibreTexts, formula 15.5.3 <https://phys.libretexts.org/Bookshelves/Electricity_and_Magnetism/Electricity_and_Magnetism_(Tatum)/15%3A_Maxwell's_Equations/15.05%3A_Maxwell's_Third_Equation>`__.

..
    TODO: shorten file name
"""

from sympy import Eq, evaluate, Expr
from symplyphysics import Quantity, validate_input, validate_output, symbols
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

magnetic_field = clone_as_vector_function(
    symbols.magnetic_field_strength,
    (position_vector, time),
)
"""
Vector of the magnetic field as a function of :attr:`~position_vector` and :attr:`~time`. See
:symbols:~magnetic_field_strength`.
"""

electric_displacement = clone_as_vector_function(
    symbols.electric_displacement,
    (position_vector, time),
)
"""
Vector of the :symbols:`electric_displacement` field as a function of :attr:`~position_vector` and
:attr:`~time`.
"""

free_current_density = clone_as_vector_function(
    symbols.current_density,
    (position_vector, time),
    display_symbol="J_f",
    display_latex="{\\vec J}_\\text{f}",
)
"""
Vector of the free (i.e. unbound) :symbols:`current_density` field as a function of
:attr:`~position_vector` and :attr:`~time`.
"""

with evaluate(False):
    _h = magnetic_field(position_vector, time)
    _curl_h = VectorCurl(_h, evaluate=False)

    _jf = free_current_density(position_vector, time)

    _d = electric_displacement(position_vector, time)
    _dd_dt = VectorDerivative(_d, time)

law = Eq(_curl_h, _jf + _dd_dt)
"""
..
    NOTE: code printers have not been implemented for `VectorCurl` yet

:code:`curl(H(r, t)) = J_f(r, t) + Derivative(D(r, t), t)`

Latex:
    .. math::
        \\text{curl} \\, \\vec H \\! \\left( \\vec r, t \\right)
            = {\\vec J}_\\text{f} \\! \\left( \\vec r, t \\right)
            + \\frac{\\partial}{\\partial t} \\vec D \\! \\left( \\vec r, t \\right)
"""


@validate_input(time_=time)
@validate_output(free_current_density)
def calculate_conductivity_current_density_at_point(
    magnetic_intensity_: Expr,
    electric_induction_: Expr,
    point_: AppliedPoint,
    time_: Quantity,
) -> QuantityCoordinateVector:
    assert_quantity_point(point_, "calculate_conductivity_current_density_at_point")

    expr = solve_for_vector(law, _jf)

    result = expr.subs({
        _h: magnetic_intensity_,
        _d: electric_induction_,
    }).doit()

    result = result.subs(point_.coordinates)
    result = result.subs(time, time_)

    result = CoordinateVector.from_expr(result)
    return QuantityCoordinateVector(result.components, result.system, result.point)

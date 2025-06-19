"""
Divergence of electric displacement field is volumetric charge density
======================================================================

The divergence of the electric induction field (also known as the electric displacement field) is
equal to the volumetric charge density at all points in space. Another form of this law states
that there exist electric charges.

**Links:**

#. `Wikipedia, first line in table <https://en.wikipedia.org/wiki/Maxwell%27s_equations#Macroscopic_formulation>`__.
#. `Physics LibreTexts, formula 15.2.3 <https://phys.libretexts.org/Bookshelves/Electricity_and_Magnetism/Electricity_and_Magnetism_(Tatum)/15%3A_Maxwell's_Equations/15.02%3A_Maxwell's_First_Equation>`__.

..
    TODO: rename file
"""

from sympy import Eq, Expr
from symplyphysics import Quantity, validate_output, symbols, clone_as_function

from symplyphysics.core.experimental.approx import assert_quantity_point
from symplyphysics.core.experimental.vectors import (clone_as_vector_function,
    clone_as_vector_symbol)
from symplyphysics.core.experimental.operators import VectorDivergence
from symplyphysics.core.experimental.coordinate_systems import AppliedPoint, CoordinateScalar

position_vector = clone_as_vector_symbol(symbols.distance_to_origin)
"""
Position vector of a point in space. See :symbols:`distance_to_origin`.
"""

electric_displacement = clone_as_vector_function(
    symbols.electric_displacement,
    (position_vector,),
)
"""
Vector field of the :symbols:`electric_displacement` as a function of the
:attr:`~position_vector`.
"""

volumetric_charge_density = clone_as_function(
    symbols.volumetric_charge_density,
    (position_vector,),
)
"""
Scalar field of the :symbols:`volumetric_charge_density` as a function of the
:attr:`~position_vector`.
"""

law = Eq(
    VectorDivergence(electric_displacement(position_vector), evaluate=False),
    volumetric_charge_density(position_vector),
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_output(volumetric_charge_density)
def calculate_charge_volumetric_density_at_point(
    electric_induction_: Expr,
    point_: AppliedPoint,
) -> Quantity:
    assert_quantity_point(point_, "calculate_charge_volumetric_density_at_point")

    result = law.lhs.subs(
        electric_displacement(position_vector),
        electric_induction_,
    ).doit()

    for base_scalar, coordinate in point_.coordinates.items():
        result = result.subs(base_scalar, coordinate)

    for scalar in result.atoms(CoordinateScalar):
        result = result.subs(scalar, scalar.scalar)

    return Quantity(result)

from typing import Optional
from pytest import approx
from sympy import Expr
from sympy.physics.units import Dimension
from sympy.physics.units.systems import SI
from symplyphysics.core.symbols.quantities import Quantity

APPROX_RELATIVE_TOLERANCE = 0.001


def approx_equal_numbers(lhs: float, rhs: float, *,
    tolerance: float = APPROX_RELATIVE_TOLERANCE) -> bool:
    rhs_approx = approx(rhs, rel=tolerance, abs=abs(lhs * tolerance))
    return lhs == rhs_approx


def approx_equal_quantities(lhs: Quantity,
    rhs: Expr | Quantity,
    *,
    tolerance: float = APPROX_RELATIVE_TOLERANCE,
    dimension: Optional[Dimension] = None) -> bool:
    rhs_quantity = rhs if isinstance(rhs, Quantity) else Quantity(rhs, dimension=dimension)
    if not SI.get_dimension_system().equivalent_dims(lhs.dimension, rhs_quantity.dimension):
        raise ValueError(
            f"Dimension of '{rhs}' is {rhs_quantity.dimension}, but it should be {lhs.dimension}")
    return approx_equal_numbers(lhs.scale_factor, rhs_quantity.scale_factor, tolerance=tolerance)


# Combined with assert for better test output
def assert_equal(lhs: Quantity,
    rhs: Expr | Quantity,
    *,
    tolerance: float = APPROX_RELATIVE_TOLERANCE,
    dimension: Optional[Dimension] = None):
    rhs_quantity = rhs if isinstance(rhs, Quantity) else Quantity(rhs, dimension=dimension)
    assert approx_equal_quantities(
        lhs, rhs, tolerance=tolerance, dimension=dimension
    ), f"Expected {lhs.scale_factor} to be equal to {rhs_quantity.scale_factor} with tolerance {tolerance}"

from typing import Optional
from pytest import approx
from sympy import N, Expr, re, im
from sympy.physics.units import Dimension
from symplyphysics.core.dimensions import assert_equivalent_dimension
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
    assert_equivalent_dimension(lhs, lhs.dimension.name, "approx_equal_quantities",
        rhs_quantity.dimension)
    if not approx_equal_numbers(
            N(im(lhs.scale_factor)), N(im(rhs_quantity.scale_factor)), tolerance=tolerance):
        return False
    return approx_equal_numbers(N(re(lhs.scale_factor)),
        N(re(rhs_quantity.scale_factor)),
        tolerance=tolerance)


# Combined with assert for better test output
def assert_equal(lhs: Expr | float,
    rhs: Expr | float,
    *,
    tolerance: float = APPROX_RELATIVE_TOLERANCE,
    dimension: Optional[Dimension] = None) -> None:
    rhs_quantity = rhs if isinstance(rhs, Quantity) else Quantity(rhs, dimension=dimension)
    # do not allow to override LHS dimension
    lhs_quantity = lhs if isinstance(lhs, Quantity) else Quantity(lhs)
    assert approx_equal_quantities(
        lhs_quantity, rhs_quantity, tolerance=tolerance, dimension=dimension
    ), f"Expected {N(lhs_quantity.scale_factor)} to be equal to {N(rhs_quantity.scale_factor)} with tolerance {tolerance}"

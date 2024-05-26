from typing import Optional
from pytest import approx
from sympy import N, Expr, re, im
from sympy.physics.units import Dimension
from symplyphysics.core.dimensions import assert_equivalent_dimension
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.core.vectors.vectors import QuantityVector

APPROX_RELATIVE_TOLERANCE = 0.001
APPROX_ABSOLUTE_TOLERANCE = 1e-20


def approx_equal_numbers(lhs: float, rhs: float, *,
    tolerance: float = APPROX_RELATIVE_TOLERANCE) -> bool:
    abs_tolerance = abs(lhs * tolerance)
    # Check for zero, as relative tolerance does not work for it
    if rhs * tolerance == rhs:
        abs_tolerance = APPROX_ABSOLUTE_TOLERANCE
    if lhs * tolerance == lhs:
        abs_tolerance = APPROX_ABSOLUTE_TOLERANCE
    rhs_approx = approx(rhs, abs=abs_tolerance)
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


def assert_equal_vectors(
    lhs: QuantityVector,
    rhs: QuantityVector,
    *,
    tolerance: float = APPROX_RELATIVE_TOLERANCE,
    dimension: Optional[Dimension] = None,
) -> None:
    for l, r in zip(lhs.components, rhs.components, strict=True):
        assert_equal(l, r, tolerance=tolerance, dimension=dimension)

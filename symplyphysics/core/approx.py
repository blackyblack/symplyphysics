from typing import Optional
from pytest import approx
from sympy import N, Expr, re, im
from sympy.physics.units import Dimension
from symplyphysics.core.dimensions import assert_equivalent_dimension
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.core.vectors.vectors import QuantityVector

APPROX_RELATIVE_TOLERANCE = 0.001


def approx_equal_numbers(lhs: float,
    rhs: float,
    *,
    relative_tolerance: Optional[float] = None,
    absolute_tolerance: Optional[float] = None) -> bool:
    if relative_tolerance is None:
        relative_tolerance = APPROX_RELATIVE_TOLERANCE
    if absolute_tolerance is None:
        absolute_tolerance = abs(lhs * relative_tolerance)
    rhs_approx = approx(rhs, rel=relative_tolerance, abs=absolute_tolerance)
    return lhs == rhs_approx


def approx_equal_quantities(lhs: Quantity,
    rhs: Expr | Quantity,
    *,
    tolerance: Optional[float] = None,
    absolute_tolerance: Optional[float] = None,
    dimension: Optional[Dimension] = None) -> bool:
    rhs_quantity = rhs if isinstance(rhs, Quantity) else Quantity(rhs, dimension=dimension)
    assert_equivalent_dimension(lhs, lhs.dimension.name, "approx_equal_quantities", rhs_quantity)
    if not approx_equal_numbers(N(im(lhs.scale_factor)),
        N(im(rhs_quantity.scale_factor)),
        relative_tolerance=tolerance,
        absolute_tolerance=absolute_tolerance):
        return False
    return approx_equal_numbers(N(re(lhs.scale_factor)),
        N(re(rhs_quantity.scale_factor)),
        relative_tolerance=tolerance,
        absolute_tolerance=absolute_tolerance)


# Combined with assert for better test output
def assert_equal(lhs: Expr | float,
    rhs: Expr | float,
    *,
    tolerance: Optional[float] = None,
    absolute_tolerance: Optional[float] = None,
    dimension: Optional[Dimension] = None) -> None:
    rhs_quantity = rhs if isinstance(rhs, Quantity) else Quantity(rhs, dimension=dimension)
    # do not allow to override LHS dimension
    lhs_quantity = lhs if isinstance(lhs, Quantity) else Quantity(lhs)
    expected_tolerance = tolerance if absolute_tolerance is None else absolute_tolerance
    expected_tolerance = APPROX_RELATIVE_TOLERANCE if expected_tolerance is None else expected_tolerance
    assert approx_equal_quantities(
        lhs_quantity,
        rhs_quantity,
        tolerance=tolerance,
        absolute_tolerance=absolute_tolerance,
        dimension=dimension
    ), f"Expected {N(lhs_quantity.scale_factor)} to be equal to {N(rhs_quantity.scale_factor)} with tolerance {expected_tolerance}"


def assert_equal_vectors(
    lhs: QuantityVector,
    rhs: QuantityVector,
    *,
    tolerance: Optional[float] = None,
    absolute_tolerance: Optional[float] = None,
    dimension: Optional[Dimension] = None,
) -> None:
    for l, r in zip(lhs.components, rhs.components, strict=True):
        assert_equal(l,
            r,
            tolerance=tolerance,
            absolute_tolerance=absolute_tolerance,
            dimension=dimension)

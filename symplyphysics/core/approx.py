"""
This module provides the functionality for checking and asserting that two numbers, quantities or
quantity vectors are equal to each other within some tolerance.

* `approx_equal_numbers` checks if two numbers are equal.
* `approx_equal_quantities` checks if two quantities are equal.
* `assert_equal` asserts that two quantities are equal.
* `assert_equal_vectors` asserts that two quantity vectors are equal.
"""

from typing import Optional
from pytest import approx
from sympy import N, re, im
from sympy.physics.units import Dimension
from symplyphysics.core.dimensions import assert_equivalent_dimension, ScalarValue
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.core.vectors.vectors import QuantityVector

_APPROX_RELATIVE_TOLERANCE = 0.001


def approx_equal_numbers(
    lhs: float,
    rhs: float,
    *,
    relative_tolerance: Optional[float] = None,
    absolute_tolerance: Optional[float] = None,
) -> bool:
    """
    Checks if the floats ``lhs`` and ``rhs`` are equal within ``relative_tolerance`` and
    ``absolute_tolerance``.

    For more information, refer to the documentation of `pytest.approx`.
    """

    if relative_tolerance is None:
        relative_tolerance = _APPROX_RELATIVE_TOLERANCE
    if absolute_tolerance is None:
        absolute_tolerance = abs(lhs * relative_tolerance)

    rhs_approx = approx(rhs, rel=relative_tolerance, abs=absolute_tolerance)
    return lhs == rhs_approx


def approx_equal_quantities(
    lhs: Quantity,
    rhs: ScalarValue,
    *,
    relative_tolerance: Optional[float] = None,
    absolute_tolerance: Optional[float] = None,
    dimension: Optional[Dimension] = None,
) -> bool:
    """
    Checks if the quantity ``lhs`` is equal to ``rhs`` within ``relative_tolerance`` and
    ``absolute_tolerance``.

    Note that ``rhs`` could be a number, in which case ``dimension`` should be set to the intended
    dimension of ``rhs``.
    """

    if not isinstance(rhs, Quantity):
        rhs = Quantity(rhs, dimension=dimension)

    assert_equivalent_dimension(lhs, lhs.dimension.name, "approx_equal_quantities", rhs)

    im_condition = approx_equal_numbers(
        float(im(lhs.scale_factor)),
        float(im(rhs.scale_factor)),
        relative_tolerance=relative_tolerance,
        absolute_tolerance=absolute_tolerance,
    )

    return im_condition and approx_equal_numbers(
        float(re(lhs.scale_factor)),
        float(re(rhs.scale_factor)),
        relative_tolerance=relative_tolerance,
        absolute_tolerance=absolute_tolerance,
    )


# Combined with assert for better test output
def assert_equal(
    lhs: ScalarValue,
    rhs: ScalarValue,
    *,
    relative_tolerance: Optional[float] = None,
    absolute_tolerance: Optional[float] = None,
    dimension: Optional[Dimension] = None,
) -> None:
    """
    Asserts that ``lhs`` and ``rhs`` are equal to each other within ``relative_tolerance`` and/or
    ``absolute_tolerance``.

    Note that ``rhs`` could be a number, in which case ``dimension`` should be set to the intended
    dimension of ``rhs``. The dimension of ``lhs`` is left untouched.
    """

    if not isinstance(rhs, Quantity):
        rhs = Quantity(rhs, dimension=dimension)

    if not isinstance(lhs, Quantity):
        # Do not allow to override LHS dimension
        lhs = Quantity(lhs)

    if absolute_tolerance is not None:
        expected_tolerance = absolute_tolerance
    elif relative_tolerance is not None:
        expected_tolerance = relative_tolerance
    else:
        expected_tolerance = _APPROX_RELATIVE_TOLERANCE

    def error_message() -> str:
        return f"Expected {N(lhs.scale_factor)} to be equal to {N(rhs.scale_factor)} with tolerance {expected_tolerance}"

    assert approx_equal_quantities(
        lhs,
        rhs,
        relative_tolerance=relative_tolerance,
        absolute_tolerance=absolute_tolerance,
        dimension=dimension,
    ), error_message()


def assert_equal_vectors(
    lhs: QuantityVector,
    rhs: QuantityVector,
    *,
    relative_tolerance: Optional[float] = None,
    absolute_tolerance: Optional[float] = None,
    dimension: Optional[Dimension] = None,
) -> None:
    """
    Asserts that the quantity vectors ``lhs`` and ``rhs`` are equal to each other component-wise
    within ``relative tolerance`` and/or ``absolute_tolerance``.
    """

    for left_component, right_component in zip(lhs.components, rhs.components, strict=True):
        assert_equal(
            left_component,
            right_component,
            relative_tolerance=relative_tolerance,
            absolute_tolerance=absolute_tolerance,
            dimension=dimension,
        )


__all__ = [
    "approx_equal_numbers",
    "approx_equal_quantities",
    "assert_equal",
    "assert_equal_vectors",
]

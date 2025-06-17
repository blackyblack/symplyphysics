from typing import Any, Optional
from sympy.physics.units import Dimension

from symplyphysics.core.approx import assert_equal
from symplyphysics.core.dimensions.dimensions import assert_equivalent_dimension

from .coordinate_systems import QuantityCoordinateVector, AppliedPoint


def assert_equal_vectors(
    lhs: Any,
    rhs: Any,
    *,
    relative_tolerance: Optional[float] = None,
    absolute_tolerance: Optional[float] = None,
    dimension: Optional[Dimension] = None,
) -> None:
    if isinstance(lhs, QuantityCoordinateVector) and isinstance(rhs, QuantityCoordinateVector):
        assert lhs.system == rhs.system
        assert lhs.point == rhs.point

        for a, b in zip(lhs.components, rhs.components):
            assert_equal(
                a,
                b,
                relative_tolerance=relative_tolerance,
                absolute_tolerance=absolute_tolerance,
                dimension=dimension,
            )

        return

    assert lhs == 0 and rhs == 0


def assert_quantity_point(point: AppliedPoint, func: Optional[str] = None) -> None:
    func = func or "assert_quantity_point"

    for base_scalar, coordinate in point.coordinates.items():
        assert_equivalent_dimension(
            coordinate,
            f"point_{base_scalar.display_name}",
            func,
            base_scalar.dimension,
        )


__all__ = [
    "assert_equal_vectors",
    "assert_quantity_point",
]

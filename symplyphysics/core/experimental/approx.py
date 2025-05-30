from typing import Any, Optional
from sympy.physics.units import Dimension

from symplyphysics.core.approx import assert_equal

from .coordinate_systems import QuantityCoordinateVector


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

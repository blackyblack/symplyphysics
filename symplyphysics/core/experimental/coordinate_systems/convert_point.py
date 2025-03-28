from typing import TypeVar

from ..points import BasePoint, PointSymbol, AppliedPoint
from .coordinate_systems import BaseCoordinateSystem
from .express_base_scalars import express_base_scalars

_P = TypeVar("_P", bound=BasePoint)


def convert_point(point: _P, new_system: BaseCoordinateSystem) -> _P:
    if isinstance(point, PointSymbol):
        return point

    if isinstance(point, AppliedPoint):
        # Point coordinates change contravariantly
        conversion = express_base_scalars(new_system, point.system)
        new_coordinates = {
            new_scalar: expr.subs(point.coordinates) for new_scalar, expr in conversion.items()
        }
        return AppliedPoint(new_coordinates, new_system)

    raise TypeError(f"Unknown point type {type(point).__name__}.")

from typing import TypeVar

from ..points import BasePoint, PointSymbol, AppliedPoint
from ..vectors import VectorExpr
from .coordinate_systems import BaseCoordinateSystem
from .express_base_scalars import express_base_scalars
from .express_base_vectors import express_base_vectors

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


def convert_vector(
    vector: VectorExpr,
    old_point: BasePoint,
    new_system: BaseCoordinateSystem,
) -> VectorExpr:
    if isinstance(old_point, PointSymbol):
        return vector

    assert isinstance(old_point, AppliedPoint)

    new_point = convert_point(old_point, new_system)

    conversion = express_base_vectors(
        old_point.system,
        new_system,
        old_args=(old_point,),
        new_args=(new_point,),
    )

    new_vector = vector.subs(conversion)

    for new_scalar, new_coordinate in new_point.coordinates.items():
        new_vector = new_vector.subs(new_scalar, new_coordinate)

    return new_vector


__all__ = [
    "convert_point",
    "convert_vector",
]

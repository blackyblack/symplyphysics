from ..points import AppliedPoint
from ..vectors import VectorExpr
from .coordinate_systems import BaseCoordinateSystem
from .express_base_scalars import express_base_scalars
from .express_base_vectors import express_base_vectors


def convert_point(point: AppliedPoint, new_system: BaseCoordinateSystem) -> AppliedPoint:
    # Point coordinates change contravariantly
    conversion = express_base_scalars(new_system, point.system)
    new_coordinates = [expr.subs(point.coordinates) for expr in conversion.values()]
    return AppliedPoint(new_coordinates, new_system)


def convert_vector(
    vector: VectorExpr,
    old_point: AppliedPoint,
    new_system: BaseCoordinateSystem,
) -> VectorExpr:
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

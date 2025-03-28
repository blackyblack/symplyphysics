from ..vectors import VectorExpr
from ..points import BasePoint, PointSymbol, AppliedPoint
from .coordinate_systems import BaseCoordinateSystem
from .convert_point import convert_point
from .express_base_vectors import express_base_vectors


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

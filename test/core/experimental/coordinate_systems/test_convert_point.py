from symplyphysics.core.experimental.coordinate_systems import (
    CartesianCoordinateSystem,
    CylindricalCoordinateSystem,
    SphericalCoordinateSystem,
)
from symplyphysics.core.experimental.coordinate_systems.convert_point import convert_point
from symplyphysics.core.experimental.points import PointSymbol


def test_convert_point_symbol() -> None:
    p = PointSymbol("P")

    cart = CartesianCoordinateSystem()
    assert convert_point(p, cart) == p

    cyl = CylindricalCoordinateSystem()
    assert convert_point(p, cyl) == p

    sph = SphericalCoordinateSystem()
    assert convert_point(p, sph) == p


# conversion of AppliedPoint is covered in `.test_scalar_conversion.py`

from symplyphysics.core.experimental.coordinate_systems import (
    CartesianCoordinateSystem,
    CylindricalCoordinateSystem,
    SphericalCoordinateSystem,
    convert_point,
)
from symplyphysics.core.experimental.points import PointSymbol


def test_convert_point_symbol() -> None:
    # Points are independent of the choice of a coordinate system, although the coordinates of the
    # points do change (which is tested in `.test_scalar_conversion.py`).

    p = PointSymbol("P")

    cart = CartesianCoordinateSystem()
    assert convert_point(p, cart) == p

    cyl = CylindricalCoordinateSystem()
    assert convert_point(p, cyl) == p

    sph = SphericalCoordinateSystem()
    assert convert_point(p, sph) == p


# conversion of AppliedPoint is covered in `.test_scalar_conversion.py`

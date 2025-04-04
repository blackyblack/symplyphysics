from pytest import raises
from sympy import pi, Expr
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.experimental.coordinate_systems import (
    CartesianCoordinateSystem,
    CylindricalCoordinateSystem,
    SphericalCoordinateSystem,
)
from symplyphysics.core.experimental.points import AppliedPoint


def test_applied_point() -> None:
    cart_sys = CartesianCoordinateSystem()
    cyl_sys = CylindricalCoordinateSystem()
    sph_sys = SphericalCoordinateSystem()

    cart_pt = AppliedPoint([0, -1, 4], cart_sys)
    assert all(isinstance(coordinate, Expr) for coordinate in cart_pt.coordinates.values())
    assert expr_equals(cart_pt[cart_sys.x], 0)
    assert expr_equals(cart_pt[cart_sys.y], -1)
    assert expr_equals(cart_pt[cart_sys.z], 4)
    assert cart_pt == cart_pt  # pylint: disable=comparison-with-itself

    # Not enough base scalars defined
    with raises(ValueError):
        AppliedPoint([4], cart_sys)
    with raises(ValueError):
        AppliedPoint([4, 5], cart_sys)

    cyl_pt = AppliedPoint([1, pi / 3, -3], cyl_sys)
    assert all(isinstance(coordinate, Expr) for coordinate in cyl_pt.coordinates.values())
    assert expr_equals(cyl_pt[cyl_sys.rho], 1)
    assert expr_equals(cyl_pt[cyl_sys.phi], pi / 3)
    assert expr_equals(cyl_pt[cyl_sys.z], -3)
    assert cyl_pt != cart_pt

    sph_pt = AppliedPoint([5, pi / 2, pi], sph_sys)
    assert all(isinstance(coordinate, Expr) for coordinate in sph_pt.coordinates.values())
    assert expr_equals(sph_pt[sph_sys.r], 5)
    assert expr_equals(sph_pt[sph_sys.theta], pi / 2)
    assert expr_equals(sph_pt[sph_sys.phi], pi)
    assert sph_pt not in (cart_pt, cyl_pt)

from pytest import raises
from sympy import pi, Expr
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.experimental.coordinate_systems import (
    CartesianCoordinateSystem,
    CylindricalCoordinateSystem,
    SphericalCoordinateSystem,
)
from symplyphysics.core.experimental.points import BasePoint, PointSymbol, AppliedPoint, PointCoordinate


def test_point_symbol() -> None:
    p = PointSymbol()
    assert p.display_name.startswith("PT")
    assert p.display_latex.startswith("P")
    assert isinstance(p, BasePoint)
    assert p == p  # pylint: disable=comparison-with-itself

    cart_sys = CartesianCoordinateSystem()
    cyl_sys = CylindricalCoordinateSystem()

    p_x: PointCoordinate = PointCoordinate(p, cart_sys.x)
    assert p_x.point == p
    assert p_x.base_scalar == cart_sys.x

    p_phi: PointCoordinate = PointCoordinate(p, cyl_sys.phi)
    assert p_phi.point == p
    assert p_phi.base_scalar == cyl_sys.phi

    q = PointSymbol(display_name="Q")
    assert q.display_name == "Q"
    assert q.display_latex == "Q"
    assert isinstance(q, BasePoint)
    assert q != p

    k = PointSymbol(display_latex="K")
    assert k.display_name.startswith("PT")
    assert k.display_latex == "K"
    assert isinstance(k, BasePoint)
    assert k not in (p, q)

    xi = PointSymbol(display_name="Xi", display_latex="\\Xi")
    assert xi.display_name == "Xi"
    assert xi.display_latex == "\\Xi"
    assert isinstance(xi, BasePoint)
    assert xi not in (p, q, k)


def test_applied_point() -> None:
    cart_sys = CartesianCoordinateSystem()
    cyl_sys = CylindricalCoordinateSystem()
    sph_sys = SphericalCoordinateSystem()

    cart_pt = AppliedPoint([0, -1, 4], cart_sys)
    assert all(isinstance(coordinate, Expr) for coordinate in cart_pt.coordinates.values())
    assert expr_equals(cart_pt[cart_sys.x], 0)
    assert expr_equals(cart_pt[cart_sys.y], -1)
    assert expr_equals(cart_pt[cart_sys.z], 4)
    assert isinstance(cart_pt, BasePoint)
    assert cart_pt == cart_pt  # pylint: disable=comparison-with-itself

    # Not enough base scalars defined
    with raises(ValueError):
        AppliedPoint([4], cart_sys)
    with raises(ValueError):
        AppliedPoint([4, 5], cart_sys)

    assert PointCoordinate(cart_pt, cart_sys.x) == 0
    assert PointCoordinate(cart_pt, cart_sys.y) == -1
    assert PointCoordinate(cart_pt, cart_sys.z) == 4
    with raises(KeyError):
        assert PointCoordinate(cart_pt, cyl_sys.phi) == 4
    with raises(KeyError):
        assert PointCoordinate(cart_pt, sph_sys.r) == 4
    # evaluate=False delays the evaluation
    bad_coordinate = PointCoordinate(cart_pt, cyl_sys.rho, evaluate=False)
    with raises(KeyError):
        bad_coordinate.doit()

    cyl_pt = AppliedPoint([1, pi / 3, -3], cyl_sys)
    assert all(isinstance(coordinate, Expr) for coordinate in cyl_pt.coordinates.values())
    assert expr_equals(cyl_pt[cyl_sys.rho], 1)
    assert expr_equals(cyl_pt[cyl_sys.phi], pi / 3)
    assert expr_equals(cyl_pt[cyl_sys.z], -3)
    assert isinstance(cyl_pt, BasePoint)
    assert cyl_pt != cart_pt

    sph_pt = AppliedPoint([5, pi / 2, pi], sph_sys)
    assert all(isinstance(coordinate, Expr) for coordinate in sph_pt.coordinates.values())
    assert expr_equals(sph_pt[sph_sys.r], 5)
    assert expr_equals(sph_pt[sph_sys.theta], pi / 2)
    assert expr_equals(sph_pt[sph_sys.phi], pi)
    assert isinstance(sph_pt, BasePoint)
    assert sph_pt not in (cart_pt, cyl_pt)

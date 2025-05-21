from dataclasses import dataclass
from pytest import raises, fixture
from sympy import pi, Expr, Symbol as SymSymbol
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.experimental.coordinate_systems import (CartesianCoordinateSystem,
    CylindricalCoordinateSystem, SphericalCoordinateSystem, AppliedPoint)
from symplyphysics.core.experimental.coordinate_systems.point import (check_point_with_system,
    GLOBAL_POINT)


@dataclass(frozen=True, kw_only=True)
class Args:
    cart_sys: CartesianCoordinateSystem
    cart_pt: AppliedPoint

    cyl_sys: CylindricalCoordinateSystem
    cyl_pt: AppliedPoint

    sph_sys: SphericalCoordinateSystem
    sph_pt: AppliedPoint


@fixture(name="test_args")
def test_args_fixture() -> Args:
    cart_sys = CartesianCoordinateSystem()
    cart_pt = AppliedPoint([0, -1, 4], cart_sys)

    cyl_sys = CylindricalCoordinateSystem()
    cyl_pt = AppliedPoint([1, pi / 3, -3], cyl_sys)

    sph_sys = SphericalCoordinateSystem()
    sph_pt = AppliedPoint([5, pi / 2, pi], sph_sys)

    return Args(
        cart_sys=cart_sys,
        cart_pt=cart_pt,
        cyl_sys=cyl_sys,
        cyl_pt=cyl_pt,
        sph_sys=sph_sys,
        sph_pt=sph_pt,
    )


def test_applied_point(test_args: Args) -> None:
    assert all(
        isinstance(coordinate, Expr) for coordinate in test_args.cart_pt.coordinates.values())
    x, y, z = test_args.cart_sys.base_scalars
    assert expr_equals(test_args.cart_pt[x], 0)
    assert expr_equals(test_args.cart_pt[y], -1)
    assert expr_equals(test_args.cart_pt[z], 4)
    assert test_args.cart_pt == test_args.cart_pt  # pylint: disable=comparison-with-itself

    # Not enough base scalars defined
    with raises(ValueError):
        AppliedPoint([4], test_args.cart_sys)
    with raises(ValueError):
        AppliedPoint([4, 5], test_args.cart_sys)

    assert all(isinstance(coordinate, Expr) for coordinate in test_args.cyl_pt.coordinates.values())
    rho, phi, z = test_args.cyl_sys.base_scalars
    assert expr_equals(test_args.cyl_pt[rho], 1)
    assert expr_equals(test_args.cyl_pt[phi], pi / 3)
    assert expr_equals(test_args.cyl_pt[z], -3)
    assert test_args.cyl_pt != test_args.cart_pt

    assert all(isinstance(coordinate, Expr) for coordinate in test_args.sph_pt.coordinates.values())
    r, theta, phi = test_args.sph_sys.base_scalars
    assert expr_equals(test_args.sph_pt[r], 5)
    assert expr_equals(test_args.sph_pt[theta], pi / 2)
    assert expr_equals(test_args.sph_pt[phi], pi)
    assert test_args.sph_pt not in (test_args.cart_pt, test_args.cyl_pt)


def test_check_point_with_system(test_args: Args) -> None:
    p = SymSymbol("P")

    assert check_point_with_system(test_args.cart_sys, None) == GLOBAL_POINT
    assert check_point_with_system(test_args.cart_sys, p) == p
    assert check_point_with_system(test_args.cart_sys, test_args.cart_pt) == test_args.cart_pt

    with raises(ValueError):
        _ = check_point_with_system(test_args.cyl_sys, None)
    assert p == check_point_with_system(test_args.cyl_sys, p)
    assert test_args.cyl_pt == check_point_with_system(test_args.cyl_sys, test_args.cyl_pt)

    with raises(ValueError):
        _ = check_point_with_system(test_args.sph_sys, None)
    assert p == check_point_with_system(test_args.sph_sys, p)
    assert test_args.sph_pt == check_point_with_system(test_args.sph_sys, test_args.sph_pt)

    with raises(ValueError):
        check_point_with_system(test_args.cart_sys, test_args.cyl_pt)
    with raises(ValueError):
        check_point_with_system(test_args.cart_sys, test_args.sph_pt)

    with raises(ValueError):
        check_point_with_system(test_args.cyl_sys, test_args.cart_pt)
    with raises(ValueError):
        check_point_with_system(test_args.cyl_sys, test_args.sph_pt)

    with raises(ValueError):
        check_point_with_system(test_args.sph_sys, test_args.cart_pt)
    with raises(ValueError):
        check_point_with_system(test_args.sph_sys, test_args.cyl_pt)

from dataclasses import dataclass
from pytest import raises, fixture
from sympy import Expr, Matrix, ImmutableMatrix, Basic
from symplyphysics.core.symbols.symbols import BasicSymbol
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.experimental.vectors import VectorSymbol
from symplyphysics.core.experimental.coordinate_systems import (CartesianCoordinateSystem,
    CylindricalCoordinateSystem, SphericalCoordinateSystem, AppliedPoint, CoordinateScalar,
    CoordinateVector)
from symplyphysics.core.experimental.coordinate_systems.point import GLOBAL_POINT


@dataclass(frozen=True, kw_only=True)
class Args:
    p: BasicSymbol

    cart_sys: CartesianCoordinateSystem
    cart_pt: AppliedPoint

    cyl_sys: CylindricalCoordinateSystem
    cyl_pt: AppliedPoint

    sph_sys: SphericalCoordinateSystem
    sph_pt: AppliedPoint


@fixture(name="test_args")
def test_args_fixture() -> Args:
    p = BasicSymbol("P")
    cart_sys = CartesianCoordinateSystem()
    cart_pt = AppliedPoint([1, 1, 1], cart_sys)

    cyl_sys = CylindricalCoordinateSystem()
    cyl_pt = AppliedPoint([1, 1, 1], cyl_sys)

    sph_sys = SphericalCoordinateSystem()
    sph_pt = AppliedPoint([1, 1, 1], sph_sys)

    return Args(
        p=p,
        cart_sys=cart_sys,
        cart_pt=cart_pt,
        cyl_sys=cyl_sys,
        cyl_pt=cyl_pt,
        sph_sys=sph_sys,
        sph_pt=sph_pt,
    )


# pylint: disable-next=too-many-statements
def test_scalar(test_args: Args) -> None:
    assert CoordinateScalar(0, test_args.cart_sys) == 0
    assert CoordinateScalar(0, test_args.cart_sys, test_args.p) == 0
    assert CoordinateScalar(0, test_args.cart_sys, test_args.cart_pt) == 0

    cart_scalar_none = CoordinateScalar(1, test_args.cart_sys)
    assert isinstance(cart_scalar_none, CoordinateScalar)
    assert expr_equals(cart_scalar_none.scalar, 1)
    assert isinstance(cart_scalar_none.scalar, Expr)
    assert cart_scalar_none.system == test_args.cart_sys  # NOTE: same point for all scalars in a Cartesian system
    assert cart_scalar_none.point == GLOBAL_POINT

    cart_scalar_p = CoordinateScalar(1, test_args.cart_sys, test_args.p)
    assert isinstance(cart_scalar_p, CoordinateScalar)
    assert expr_equals(cart_scalar_p.scalar, 1)
    assert isinstance(cart_scalar_p.scalar, Expr)
    assert cart_scalar_p.system == test_args.cart_sys
    assert cart_scalar_p.point == test_args.p
    assert cart_scalar_p != 1

    assert cart_scalar_p != cart_scalar_none

    cart_scalar_q = CoordinateScalar(1, test_args.cart_sys, test_args.cart_pt)
    assert isinstance(cart_scalar_q, CoordinateScalar)
    assert expr_equals(cart_scalar_q.scalar, 1)
    assert isinstance(cart_scalar_q.scalar, Expr)
    assert cart_scalar_q.system == test_args.cart_sys
    assert cart_scalar_q.point == test_args.cart_pt
    assert cart_scalar_q != 1

    assert cart_scalar_q != cart_scalar_none
    assert cart_scalar_q != cart_scalar_p

    with raises(ValueError):
        _ = CoordinateScalar(1, test_args.cyl_sys)

    cyl_scalar_p = CoordinateScalar(1, test_args.cyl_sys, test_args.p)
    assert isinstance(cyl_scalar_p, CoordinateScalar)
    assert expr_equals(cyl_scalar_p.scalar, 1)
    assert isinstance(cyl_scalar_p.scalar, Expr)
    assert cyl_scalar_p.system == test_args.cyl_sys
    assert cyl_scalar_p.point == test_args.p

    assert cyl_scalar_p != 1
    assert cyl_scalar_p != cart_scalar_none
    assert cyl_scalar_p != cart_scalar_p
    assert cyl_scalar_p != cart_scalar_q

    cyl_scalar_q = CoordinateScalar(1, test_args.cyl_sys, test_args.cyl_pt)
    assert isinstance(cyl_scalar_q, CoordinateScalar)
    assert expr_equals(cyl_scalar_q.scalar, 1)
    assert isinstance(cyl_scalar_q.scalar, Expr)
    assert cyl_scalar_q.system == test_args.cyl_sys
    assert cyl_scalar_q.point == test_args.cyl_pt

    assert cyl_scalar_q != 1
    assert cyl_scalar_q != cyl_scalar_p
    assert cyl_scalar_q != cart_scalar_none
    assert cyl_scalar_q != cart_scalar_p
    assert cyl_scalar_q != cart_scalar_q

    with raises(ValueError):
        _ = CoordinateScalar(1, test_args.sph_sys)

    sph_scalar_p = CoordinateScalar(1, test_args.sph_sys, test_args.p)
    assert isinstance(sph_scalar_p, CoordinateScalar)
    assert expr_equals(sph_scalar_p.scalar, 1)
    assert isinstance(sph_scalar_p.scalar, Expr)
    assert sph_scalar_p.system == test_args.sph_sys
    assert sph_scalar_p.point == test_args.p

    assert sph_scalar_p != 1
    assert sph_scalar_p != cyl_scalar_p

    sph_scalar_q = CoordinateScalar(1, test_args.sph_sys, test_args.sph_pt)
    assert isinstance(sph_scalar_q, CoordinateScalar)
    assert expr_equals(sph_scalar_q.scalar, 1)
    assert isinstance(sph_scalar_q.scalar, Expr)
    assert sph_scalar_q.system == test_args.sph_sys
    assert sph_scalar_q.point == test_args.sph_pt

    assert sph_scalar_q != 1
    assert sph_scalar_q != sph_scalar_p


def test_scalar_bad_input(test_args: Args) -> None:
    with raises(ValueError):
        CoordinateScalar("1", test_args.cart_sys)

    with raises(TypeError):
        CoordinateScalar(Basic(), test_args.cart_sys)
    with raises(TypeError):
        CoordinateScalar(int, test_args.cart_sys)

    # matrix input
    with raises(ValueError):
        CoordinateScalar(Matrix([1, 1, 1]), test_args.cart_sys)
    with raises(ValueError):
        CoordinateScalar(ImmutableMatrix([1, 1, 1]), test_args.cart_sys)

    # vector input
    v = VectorSymbol("v")
    w = VectorSymbol("w")
    with raises(ValueError):
        CoordinateScalar(v, test_args.cart_sys)
    with raises(ValueError):
        CoordinateScalar(v + w, test_args.cart_sys)
    with raises(ValueError):
        CoordinateScalar(CoordinateVector([1, 1, 1], test_args.cart_sys), test_args.cart_sys)

    # point and system mismatch
    with raises(ValueError):
        CoordinateScalar(1, test_args.cart_sys, test_args.cyl_pt)
    with raises(ValueError):
        CoordinateScalar(1, test_args.cart_sys, test_args.sph_pt)
    with raises(ValueError):
        CoordinateScalar(1, test_args.cyl_sys, test_args.cart_pt)
    with raises(ValueError):
        CoordinateScalar(1, test_args.cyl_sys, test_args.sph_pt)
    with raises(ValueError):
        CoordinateScalar(1, test_args.sph_sys, test_args.cart_pt)
    with raises(ValueError):
        CoordinateScalar(1, test_args.sph_sys, test_args.cyl_pt)

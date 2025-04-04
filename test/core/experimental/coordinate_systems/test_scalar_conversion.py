from dataclasses import dataclass
from pytest import fixture, raises
from sympy import S, pi, sqrt
from symplyphysics.core.experimental.points import AppliedPoint
from symplyphysics.core.experimental.coordinate_systems import (
    CartesianCoordinateSystem,
    CylindricalCoordinateSystem,
    SphericalCoordinateSystem,
    express_base_scalars,
    convert_point,
)


@dataclass
class Args:
    cart: CartesianCoordinateSystem
    cyl: CylindricalCoordinateSystem
    sph: SphericalCoordinateSystem

    p_cart: AppliedPoint
    p_cyl: AppliedPoint
    p_sph: AppliedPoint


@fixture(name="test_args")
def test_args_fixture() -> Args:
    cart = CartesianCoordinateSystem()
    cyl = CylindricalCoordinateSystem()
    sph = SphericalCoordinateSystem()

    # These represent the same physical point in different coordinate systems
    p_cart = AppliedPoint({cart.x: S(-1), cart.y: S(0), cart.z: S(-1)}, cart)
    p_cyl = AppliedPoint({cyl.rho: S(1), cyl.phi: pi, cyl.z: S(-1)}, cyl)
    p_sph = AppliedPoint({sph.r: sqrt(2), sph.theta: 3 * pi / 4, sph.phi: pi}, sph)

    return Args(cart=cart, cyl=cyl, sph=sph, p_cart=p_cart, p_cyl=p_cyl, p_sph=p_sph)


equal_points = AppliedPoint.equals


def test_cartesian_to_cartesian(test_args: Args) -> None:
    new_cart = CartesianCoordinateSystem()

    new_point = convert_point(test_args.p_cart, new_cart)

    correct_point = AppliedPoint.from_iterable(test_args.p_cart.coordinates.values(), new_cart)
    assert equal_points(new_point, correct_point)

    wrong_point = AppliedPoint.from_iterable([pi / 2] * 3, new_cart)
    assert not equal_points(new_point, wrong_point)


def test_cartesian_to_cylindrical(test_args: Args) -> None:
    new_point = convert_point(test_args.p_cart, test_args.cyl)

    assert equal_points(new_point, test_args.p_cyl)

    wrong_point = AppliedPoint.from_iterable([pi / 2] * 3, test_args.cyl)
    assert not equal_points(new_point, wrong_point)


def test_cartesian_to_spherical(test_args: Args) -> None:
    new_point = convert_point(test_args.p_cart, test_args.sph)

    assert equal_points(new_point, test_args.p_sph)

    wrong_point = AppliedPoint.from_iterable([pi / 2] * 3, test_args.sph)
    assert not equal_points(new_point, wrong_point)


def test_cylindrical_to_cartesian(test_args: Args) -> None:
    new_point = convert_point(test_args.p_cyl, test_args.cart)

    assert equal_points(new_point, test_args.p_cart)

    wrong_point = AppliedPoint.from_iterable([pi / 2] * 3, test_args.cart)
    assert not equal_points(new_point, wrong_point)


def test_cylindrical_to_cylindrical(test_args: Args) -> None:
    new_cyl = CylindricalCoordinateSystem()

    new_point = convert_point(test_args.p_cyl, new_cyl)

    correct_point = AppliedPoint.from_iterable(test_args.p_cyl.coordinates.values(), new_cyl)
    assert equal_points(new_point, correct_point)

    wrong_point = AppliedPoint.from_iterable([pi / 2] * 3, new_cyl)
    assert not equal_points(new_point, wrong_point)


def test_cylindrical_to_spherical(test_args: Args) -> None:
    new_point = convert_point(test_args.p_cyl, test_args.sph)

    assert equal_points(new_point, test_args.p_sph)

    wrong_point = AppliedPoint.from_iterable([pi / 2] * 3, test_args.sph)
    assert not equal_points(new_point, wrong_point)


def test_spherical_to_cartesian(test_args: Args) -> None:
    new_point = convert_point(test_args.p_sph, test_args.cart)

    assert equal_points(new_point, test_args.p_cart)

    wrong_point = AppliedPoint.from_iterable([pi / 2] * 3, test_args.cart)
    assert not equal_points(new_point, wrong_point)


def test_spherical_to_cylindrical(test_args: Args) -> None:
    new_point = convert_point(test_args.p_sph, test_args.cyl)

    assert equal_points(new_point, test_args.p_cyl)

    wrong_point = AppliedPoint.from_iterable([pi / 2] * 3, test_args.cyl)
    assert not equal_points(new_point, wrong_point)


def test_spherical_to_spherical(test_args: Args) -> None:
    new_sph = SphericalCoordinateSystem()

    new_point = convert_point(test_args.p_sph, new_sph)

    correct_point = AppliedPoint.from_iterable(test_args.p_sph.coordinates.values(), new_sph)
    assert equal_points(new_point, correct_point)

    wrong_point = AppliedPoint.from_iterable([pi / 2] * 3, new_sph)
    assert not equal_points(new_point, wrong_point)


def test_unregistered_coordinate_system(test_args: Args) -> None:

    class NewCartesianCoordinateSystem(CartesianCoordinateSystem):
        pass

    new_cart = NewCartesianCoordinateSystem()

    with raises(TypeError):
        express_base_scalars(new_cart, test_args.cart)

    with raises(TypeError):
        express_base_scalars(test_args.cart, new_cart)

    # NOTE: conversion between identical coordinate systems is allowed, even if unregistered
    _ = express_base_scalars(new_cart, new_cart)


def test_same_system(test_args: Args) -> None:
    new_point = convert_point(test_args.p_cart, test_args.cart)
    assert equal_points(test_args.p_cart, new_point)

    new_point = convert_point(test_args.p_cyl, test_args.cyl)
    assert equal_points(test_args.p_cyl, new_point)

    new_point = convert_point(test_args.p_sph, test_args.sph)
    assert equal_points(test_args.p_sph, new_point)

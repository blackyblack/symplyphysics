from typing import TypeAlias
from dataclasses import dataclass
from pytest import fixture, raises
from sympy import Expr, Symbol as SymSymbol, S, pi, sqrt
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.experimental.coordinate_systems import (
    CartesianCoordinateSystem,
    CylindricalCoordinateSystem,
    SphericalCoordinateSystem,
)
from symplyphysics.core.experimental.coordinate_systems.convert_base_scalars import convert_base_scalars, ScalarConversion

Point: TypeAlias = dict[SymSymbol, Expr]


@dataclass
class Args:
    cart: CartesianCoordinateSystem
    cyl: CylindricalCoordinateSystem
    sph: SphericalCoordinateSystem

    p_cart: Point
    p_cyl: Point
    p_sph: Point


@fixture(name="test_args")
def test_args_fixture() -> Args:
    cart = CartesianCoordinateSystem()
    cyl = CylindricalCoordinateSystem()
    sph = SphericalCoordinateSystem()

    # These represent the same physical point in different coordinate systems
    p_cart = {cart.x: S(-1), cart.y: S(0), cart.z: S(-1)}
    p_cyl = {cyl.rho: S(1), cyl.phi: pi, cyl.z: S(-1)}
    p_sph = {sph.r: sqrt(2), sph.theta: 3 * pi / 4, sph.phi: pi}

    return Args(cart=cart, cyl=cyl, sph=sph, p_cart=p_cart, p_cyl=p_cyl, p_sph=p_sph)


def apply_conversion(point: Point, conversion: ScalarConversion) -> Point:
    return {new_coordinate: expr.subs(point) for (new_coordinate, expr) in conversion.items()}


def equal_points(point1: Point, point2: Point) -> bool:
    return all(expr_equals(coordinate1, point2[s]) for (s, coordinate1) in point1.items())


def test_cartesian_to_cartesian(test_args: Args) -> None:
    old_point = test_args.p_cart

    new_cart = CartesianCoordinateSystem()

    conversion = convert_base_scalars(new_cart, test_args.cart)
    new_point = apply_conversion(old_point, conversion)

    correct_point = dict(zip(new_cart.base_scalars, old_point.values()))
    assert equal_points(new_point, correct_point)

    wrong_point = dict(zip(new_cart.base_scalars, [S.Zero] * 3))
    assert not equal_points(new_point, wrong_point)


def test_cartesian_to_cylindrical(test_args: Args) -> None:
    old_point = test_args.p_cart
    conversion = convert_base_scalars(test_args.cyl, test_args.cart)
    new_point = apply_conversion(old_point, conversion)

    assert equal_points(new_point, test_args.p_cyl)

    wrong_point = dict(zip(test_args.cyl.base_scalars, [S.Zero] * 3))
    assert not equal_points(new_point, wrong_point)


def test_cartesian_to_spherical(test_args: Args) -> None:
    old_point = test_args.p_cart
    conversion = convert_base_scalars(test_args.sph, test_args.cart)
    new_point = apply_conversion(old_point, conversion)

    assert equal_points(new_point, test_args.p_sph)

    wrong_point = dict(zip(test_args.sph.base_scalars, [S.Zero] * 3))
    assert not equal_points(new_point, wrong_point)


def test_cylindrical_to_cartesian(test_args: Args) -> None:
    old_point = test_args.p_cyl
    conversion = convert_base_scalars(test_args.cart, test_args.cyl)
    new_point = apply_conversion(old_point, conversion)

    assert equal_points(new_point, test_args.p_cart)

    wrong_point = dict(zip(test_args.cart.base_scalars, [S.Zero] * 3))
    assert not equal_points(new_point, wrong_point)


def test_cylindrical_to_cylindrical(test_args: Args) -> None:
    old_point = test_args.p_cyl

    new_cyl = CylindricalCoordinateSystem()

    conversion = convert_base_scalars(new_cyl, test_args.cyl)
    new_point = apply_conversion(old_point, conversion)

    correct_point = dict(zip(new_cyl.base_scalars, old_point.values()))
    assert equal_points(new_point, correct_point)

    wrong_point = dict(zip(new_cyl.base_scalars, [S.Zero] * 3))
    assert not equal_points(new_point, wrong_point)


def test_cylindrical_to_spherical(test_args: Args) -> None:
    old_point = test_args.p_cyl
    conversion = convert_base_scalars(test_args.sph, test_args.cyl)
    new_point = apply_conversion(old_point, conversion)

    assert equal_points(new_point, test_args.p_sph)

    wrong_point = dict(zip(test_args.sph.base_scalars, [S.Zero] * 3))
    assert not equal_points(new_point, wrong_point)


def test_spherical_to_cartesian(test_args: Args) -> None:
    old_point = test_args.p_sph
    conversion = convert_base_scalars(test_args.cart, test_args.sph)
    new_point = apply_conversion(old_point, conversion)

    assert equal_points(new_point, test_args.p_cart)

    wrong_point = dict(zip(test_args.cart.base_scalars, [S.Zero] * 3))
    assert not equal_points(new_point, wrong_point)


def test_spherical_to_cylindrical(test_args: Args) -> None:
    old_point = test_args.p_sph
    conversion = convert_base_scalars(test_args.cyl, test_args.sph)
    new_point = apply_conversion(old_point, conversion)

    assert equal_points(new_point, test_args.p_cyl)

    wrong_point = dict(zip(test_args.cyl.base_scalars, [S.Zero] * 3))
    assert not equal_points(new_point, wrong_point)


def test_spherical_to_spherical(test_args: Args) -> None:
    old_point = test_args.p_sph

    new_sph = SphericalCoordinateSystem()

    conversion = convert_base_scalars(new_sph, test_args.sph)
    new_point = apply_conversion(old_point, conversion)

    correct_point = dict(zip(new_sph.base_scalars, old_point.values()))
    assert equal_points(new_point, correct_point)

    wrong_point = dict(zip(new_sph.base_scalars, [S.Zero] * 3))
    assert not equal_points(new_point, wrong_point)


def test_unregistered_coordinate_system(test_args: Args) -> None:

    class NewCartesianCoordinateSystem(CartesianCoordinateSystem):
        pass

    new_cart = NewCartesianCoordinateSystem()

    with raises(TypeError):
        convert_base_scalars(new_cart, test_args.cart)

    with raises(TypeError):
        convert_base_scalars(test_args.cart, new_cart)

    # NOTE: conversion between identical coordinate systems, even if unregistered
    _ = convert_base_scalars(new_cart, new_cart)


def test_same_system(test_args: Args) -> None:
    conversion = convert_base_scalars(test_args.cart, test_args.cart)
    old_point = test_args.p_cart
    assert equal_points(old_point, apply_conversion(old_point, conversion))

    conversion = convert_base_scalars(test_args.cyl, test_args.cyl)
    old_point = test_args.p_cyl
    assert equal_points(old_point, apply_conversion(old_point, conversion))

    conversion = convert_base_scalars(test_args.sph, test_args.sph)
    old_point = test_args.p_sph
    assert equal_points(old_point, apply_conversion(old_point, conversion))

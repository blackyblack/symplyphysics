from typing import TypeAlias
from dataclasses import dataclass
from pytest import fixture, raises
from sympy import Expr, Symbol as SymSymbol, S, pi, sqrt
from symplyphysics.core.experimental.vectors import VectorExpr, norm, ZERO
from symplyphysics.core.experimental.coordinate_systems import (
    BaseCoordinateSystem,
    CartesianCoordinateSystem,
    CylindricalCoordinateSystem,
    SphericalCoordinateSystem,
)
from symplyphysics.core.experimental.coordinate_systems.express_base_scalars import express_base_scalars
from symplyphysics.core.experimental.coordinate_systems.express_base_vectors import express_base_vectors

Point: TypeAlias = dict[SymSymbol, Expr]


@dataclass
class Args:  # pylint: disable=too-many-instance-attributes
    cart: CartesianCoordinateSystem
    cyl: CylindricalCoordinateSystem
    sph: SphericalCoordinateSystem

    p_cart: Point
    p_cyl: Point
    p_sph: Point

    v_cart: VectorExpr
    v_cyl: VectorExpr
    v_sph: VectorExpr


@fixture(name="test_args")
def test_args_fixture() -> Args:
    cart = CartesianCoordinateSystem()
    cyl = CylindricalCoordinateSystem()
    sph = SphericalCoordinateSystem()

    # These represent the same physical point in different coordinate systems
    p_cart = {cart.x: S(-1), cart.y: S(0), cart.z: S(-1)}
    p_cyl = {cyl.rho: S(1), cyl.phi: pi, cyl.z: S(-1)}
    p_sph = {sph.r: sqrt(2), sph.theta: 3 * pi / 4, sph.phi: pi}

    # The following vectors represent the same vector entity in different coordinate systems

    v_cart = cart.j + cart.k * 2

    _, e_phi_cyl, e_z = cyl.base_vectors
    v_cyl = -e_phi_cyl + e_z * 2

    e_r, e_theta, e_phi_sph = sph.base_vectors
    v_sph = e_r * -sqrt(2) - e_phi_sph + e_theta * -sqrt(2)

    return Args(
        cart=cart,
        cyl=cyl,
        sph=sph,
        p_cart=p_cart,
        p_cyl=p_cyl,
        p_sph=p_sph,
        v_cart=v_cart,
        v_cyl=v_cyl,
        v_sph=v_sph,
    )


def convert_point(
    old_point: Point,
    old_system: BaseCoordinateSystem,
    new_system: BaseCoordinateSystem,
) -> Point:
    # Point coordinates change contravariantly
    conversion = express_base_scalars(new_system, old_system)  # pylint: disable=arguments-out-of-order
    return {new_coordinate: expr.subs(old_point) for new_coordinate, expr in conversion.items()}


def convert_vector(
    old_vector: VectorExpr,
    old_system: BaseCoordinateSystem,
    new_system: BaseCoordinateSystem,
    new_point: Point,
) -> VectorExpr:
    # Basic vectors change covariantly
    conversion = express_base_vectors(old_system, new_system)
    return old_vector.subs(conversion).subs(new_point).doit()


def equal_vectors(vector1: VectorExpr, vector2: VectorExpr) -> bool:
    return bool(norm(vector1 - vector2) == 0)


def test_cartesian_to_cartesian(test_args: Args) -> None:
    new_cart = CartesianCoordinateSystem()

    new_point = convert_point(test_args.p_cart, test_args.cart, new_cart)
    new_vector = convert_vector(test_args.v_cart, test_args.cart, new_cart, new_point)

    substitution = dict(zip(test_args.cart.base_vectors, new_cart.base_vectors))
    correct_vector = test_args.v_cart.subs(substitution)

    assert equal_vectors(new_vector, correct_vector)

    assert not equal_vectors(new_vector, ZERO)


def test_cartesian_to_cylindrical(test_args: Args) -> None:
    new_vector = convert_vector(test_args.v_cart, test_args.cart, test_args.cyl, test_args.p_cyl)

    assert equal_vectors(new_vector, test_args.v_cyl)

    assert not equal_vectors(new_vector, ZERO)


def test_cartesian_to_spherical(test_args: Args) -> None:
    new_vector = convert_vector(test_args.v_sph, test_args.cart, test_args.sph, test_args.p_sph)

    assert equal_vectors(new_vector, test_args.v_sph)

    assert not equal_vectors(new_vector, ZERO)


def test_cylindrical_to_cartesian(test_args: Args) -> None:
    new_vector = convert_vector(test_args.v_cyl, test_args.cyl, test_args.cart, test_args.p_cart)

    assert equal_vectors(new_vector, test_args.v_cart)

    assert not equal_vectors(new_vector, ZERO)


def test_cylindrical_to_cylindrical(test_args: Args) -> None:
    new_cyl = CylindricalCoordinateSystem()

    new_point = convert_point(test_args.p_cyl, test_args.cyl, new_cyl)
    new_vector = convert_vector(test_args.v_cyl, test_args.cyl, new_cyl, new_point)

    substitution = dict(zip(test_args.cyl.base_vectors, new_cyl.base_vectors))
    correct_vector = test_args.v_cyl.subs(substitution)

    assert equal_vectors(new_vector, correct_vector)

    assert not equal_vectors(new_vector, ZERO)


def test_cylindrical_to_spherical(test_args: Args) -> None:
    new_vector = convert_vector(test_args.v_cyl, test_args.cyl, test_args.sph, test_args.p_sph)

    assert equal_vectors(new_vector, test_args.v_sph)

    assert not equal_vectors(new_vector, ZERO)


def test_spherical_to_cartesian(test_args: Args) -> None:
    new_vector = convert_vector(test_args.v_sph, test_args.sph, test_args.cart, test_args.p_cart)

    assert equal_vectors(new_vector, test_args.v_cart)

    assert not equal_vectors(new_vector, ZERO)


def test_spherical_to_cylindrical(test_args: Args) -> None:
    new_vector = convert_vector(test_args.v_sph, test_args.sph, test_args.cyl, test_args.p_cyl)

    assert equal_vectors(new_vector, test_args.v_cyl)

    assert not equal_vectors(new_vector, ZERO)


def test_spherical_to_spherical(test_args: Args) -> None:
    new_sph = SphericalCoordinateSystem()

    new_point = convert_point(test_args.p_sph, test_args.sph, new_sph)
    new_vector = convert_vector(test_args.v_sph, test_args.sph, new_sph, new_point)

    substitution = dict(zip(test_args.sph.base_vectors, new_sph.base_vectors))
    correct_vector = test_args.v_sph.subs(substitution)

    assert equal_vectors(new_vector, correct_vector)

    assert not equal_vectors(new_vector, ZERO)


def test_unregistered_coordinate_system(test_args: Args) -> None:

    class NewCartesianCoordinateSystem(CartesianCoordinateSystem):
        pass

    new_cart = NewCartesianCoordinateSystem()

    with raises(TypeError):
        express_base_vectors(new_cart, test_args.cart)

    with raises(TypeError):
        express_base_vectors(test_args.cart, new_cart)

    # NOTE: conversion between identical coordinate systems is allowed, even if unregistered
    _ = express_base_vectors(new_cart, new_cart)
    _ = express_base_vectors(new_cart, NewCartesianCoordinateSystem())
    _ = express_base_vectors(NewCartesianCoordinateSystem(), new_cart)


def test_same_system(test_args: Args) -> None:
    new_vector = convert_vector(test_args.v_cart, test_args.cart, test_args.cart, test_args.p_cart)
    assert equal_vectors(new_vector, test_args.v_cart)

    new_vector = convert_vector(test_args.v_cyl, test_args.cyl, test_args.cyl, test_args.p_cyl)
    assert equal_vectors(new_vector, test_args.v_cyl)

    new_vector = convert_vector(test_args.v_sph, test_args.sph, test_args.sph, test_args.p_sph)
    assert equal_vectors(new_vector, test_args.v_sph)

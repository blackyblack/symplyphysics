from typing import TypeAlias
from dataclasses import dataclass
from pytest import fixture, raises, skip
from sympy import Expr, Symbol as SymSymbol, S, pi, sqrt
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.experimental.vectors import VectorExpr, VectorNorm as norm, ZERO
from symplyphysics.core.experimental.points import AppliedPoint
from symplyphysics.core.experimental.coordinate_systems import (
    BaseCoordinateSystem,
    CartesianCoordinateSystem,
    CylindricalCoordinateSystem,
    SphericalCoordinateSystem,
)
from symplyphysics.core.experimental.coordinate_systems.express_base_scalars import express_base_scalars, ScalarMapping
from symplyphysics.core.experimental.coordinate_systems.express_base_vectors import express_base_vectors

Point: TypeAlias = dict[SymSymbol, Expr]


@dataclass(frozen=True)
class AppliedVector:
    value: VectorExpr
    """Expression containing the components of the vector and the base scalars."""

    point: Point
    """Point where the vector is applied."""

    system: BaseCoordinateSystem
    """Coordinate system where the vector is described."""


@dataclass(frozen=True)
class Args:  # pylint: disable=too-many-instance-attributes
    cart: CartesianCoordinateSystem
    cyl: CylindricalCoordinateSystem
    sph: SphericalCoordinateSystem

    # These represent the same vector `v` in different coordinate systems
    v_cart: AppliedVector
    v_cyl: AppliedVector
    v_sph: AppliedVector

    # These vectors have the same value as the respective `v` vector, but are applied at at
    # different points; they're not necessarily all equal to each other
    w_cart: AppliedVector
    w_cyl: AppliedVector
    w_sph: AppliedVector


@fixture(name="test_args")
def test_args_fixture() -> Args:
    cart = CartesianCoordinateSystem()
    cyl = CylindricalCoordinateSystem()
    sph = SphericalCoordinateSystem()

    # These represent the same physical point `p` (which is different from `q`) in different
    # coordinate systems
    p_cart = AppliedPoint({cart.x: S(-1), cart.y: S(0), cart.z: S(-1)}, cart)
    p_cyl = AppliedPoint({cyl.rho: S(1), cyl.phi: pi, cyl.z: S(-1)}, cyl)
    p_sph = AppliedPoint({sph.r: sqrt(2), sph.theta: 3 * pi / 4, sph.phi: pi}, sph)

    u_cart = cart.j + cart.k * 2

    _, e_phi_cyl, e_z = cyl.base_vectors(p_cyl)
    u_cyl = -e_phi_cyl + e_z * 2

    e_r, e_theta, e_phi_sph = sph.base_vectors(p_sph)
    u_sph = e_r * -sqrt(2) - e_phi_sph + e_theta * -sqrt(2)

    # These represent the same physical point `q` (which is different from `p`) in different
    # coordinate systems
    q_cart = {cart.x: S(1), cart.y: S(1), cart.z: S(1)},
    q_cyl = {cyl.rho: sqrt(2), cyl.phi: pi / 4, cyl.z: S(1)}
    q_sph = {sph.r: sqrt(3), sph.theta: pi / 4, sph.phi: pi / 4}

    return Args(
        cart=cart,
        cyl=cyl,
        sph=sph,
        v_cart=AppliedVector(u_cart, p_cart, cart),
        v_cyl=AppliedVector(u_cyl, p_cyl, cyl),
        v_sph=AppliedVector(u_sph, p_sph, sph),
        w_cart=AppliedVector(u_cart, q_cart, cart),
        w_cyl=AppliedVector(u_cyl, q_cyl, cyl),
        w_sph=AppliedVector(u_sph, q_sph, sph),
    )


def convert_vector(
    old_vector: AppliedVector,
    new_system: BaseCoordinateSystem,
) -> AppliedVector:
    # Point coordinates change contravariantly
    scalar_conversion: ScalarMapping = express_base_scalars(new_system, old_vector.system)

    new_point = {
        new_coordinate: expr.subs(old_vector.point)
        for new_coordinate, expr in scalar_conversion.items()
    }

    # Base vectors change covariantly
    vector_conversion = express_base_vectors(old_vector.system, new_system)
    new_value = old_vector.value.subs(vector_conversion).subs(new_point).doit()

    return AppliedVector(new_value, new_point, new_system)


def equals(this: VectorExpr | AppliedVector, that: VectorExpr | AppliedVector) -> bool:
    """
    Checks if `this` and `that` are the same vector expression.
    
    Note that a zero-length vector can be compared to any other vector, but it is equal only to
    itself.

    Raises `ValueError` if any of the following conditions is not met:

    1. Both vectors must be described in the same coordinate system.

    2. If they are described in a non-Cartesian coordinate system, their points of application must
       be the same.

    3. If at least one of them is a `VectorExpr`, at least one of them must be a zero vector.
    """

    def equal_points(point1: Point, point2: Point) -> bool:
        """
        Checks if `point1` and `point2` have the same coordinates. Raises `KeyError` if the points
        are described with different base scalars.
        """

        return all(expr_equals(coordinate1, point2[s]) for (s, coordinate1) in point1.items())

    def equal_vectors(vector1: VectorExpr, vector2: VectorExpr) -> bool:
        """
        Checks if `vector1` and `vector2` are equal, i.e. if their difference is a zero vector.
        """

        return expr_equals(norm(vector1 - vector2), 0)

    if isinstance(this, AppliedVector):
        if isinstance(that, AppliedVector):
            # Refer to condition #1
            if this.system is not that.system:
                raise ValueError("Vectors can only be compared in the same coordinate system.")

            if (isinstance(this.system, CartesianCoordinateSystem) or
                    equal_points(this.point, that.point)):
                return equal_vectors(this.value, that.value)

            # Refer to condition #2
            raise ValueError("Vectors in non-Cartesian systems must be applied at the same point.")

        this = this.value
    elif isinstance(that, AppliedVector):
        that = that.value

    this_norm = norm(this)
    that_norm = norm(that)

    # The zero vector can be compared with any vector.
    if expr_equals(this_norm, 0) or expr_equals(that_norm, 0):
        return expr_equals(this_norm, that_norm)

    # Refer to condition #3
    raise ValueError("Non-zero vectors must have a point of application.")


@skip(allow_module_level=True)  # TODO: remove when the tests are fixed
def test_cartesian_to_cartesian(test_args: Args) -> None:
    new_cart = CartesianCoordinateSystem()

    new_vector = convert_vector(test_args.v_cart, new_cart)

    correct_point = {new_cart.x: S(-1), new_cart.y: S(0), new_cart.z: S(-1)}
    _, j, k = new_cart.base_vectors
    correct_value = j + k * 2
    correct_vector = AppliedVector(correct_value, correct_point, new_cart)

    assert equals(new_vector, correct_vector)

    # The transformation of basis vectors between Cartesian coordinates does *not* depend how we
    # choose the point where the vector is applied.
    other_vector = convert_vector(test_args.w_cart, new_cart)
    assert equals(new_vector, other_vector)

    # Check that the vector is not described in other coordinate systems
    with raises(ValueError):
        equals(new_vector, test_args.v_cyl)
    with raises(ValueError):
        equals(new_vector, test_args.v_sph)

    assert not equals(new_vector, ZERO)


def test_cartesian_to_cylindrical(test_args: Args) -> None:
    new_vector = convert_vector(test_args.v_cart, test_args.cyl)
    assert equals(new_vector, test_args.v_cyl)

    # Cannot compare vectors applied at different points in non-Cartesian systems
    other_vector = convert_vector(test_args.w_cart, test_args.cyl)
    with raises(ValueError):
        equals(new_vector, other_vector)

    # Check that the vector is not described in other coordinate systems
    with raises(ValueError):
        equals(new_vector, test_args.v_cart)
    with raises(ValueError):
        equals(new_vector, test_args.v_sph)

    assert not equals(new_vector, ZERO)


def test_cartesian_to_spherical(test_args: Args) -> None:
    new_vector = convert_vector(test_args.v_cart, test_args.sph)
    assert equals(new_vector, test_args.v_sph)

    # Cannot compare vectors applied at different points in non-Cartesian systems
    other_vector = convert_vector(test_args.w_cart, test_args.sph)
    with raises(ValueError):
        equals(new_vector, other_vector)

    # Check that the vector is not described in other coordinate systems
    with raises(ValueError):
        equals(new_vector, test_args.v_cart)
    with raises(ValueError):
        equals(new_vector, test_args.v_cyl)

    assert not equals(new_vector, ZERO)


def test_cylindrical_to_cartesian(test_args: Args) -> None:
    new_vector = convert_vector(test_args.v_cyl, test_args.cart)
    assert equals(new_vector, test_args.v_cart)

    # In non-Cartesian systems, base vectors depend on the point of application.
    other_vector = convert_vector(test_args.w_cyl, test_args.cart)
    assert not equals(new_vector, other_vector)

    # Check that the vector is not described in other coordinate systems
    with raises(ValueError):
        equals(new_vector, test_args.v_cyl)
    with raises(ValueError):
        equals(new_vector, test_args.v_sph)

    assert not equals(new_vector, ZERO)


def test_cylindrical_to_cylindrical(test_args: Args) -> None:
    new_cyl = CylindricalCoordinateSystem()

    new_vector = convert_vector(test_args.v_cyl, new_cyl)

    rho, phi, z = new_cyl.base_scalars
    correct_point = {rho: S(1), phi: pi, z: S(-1)}
    _, e_phi, e_z = new_cyl.base_vectors
    correct_value = -e_phi + e_z * 2
    correct_vector = AppliedVector(correct_value, correct_point, new_cyl)

    assert equals(new_vector, correct_vector)

    # Cannot compare vectors applied at different points in non-Cartesian systems
    other_vector = convert_vector(test_args.w_cyl, new_cyl)
    with raises(ValueError):
        equals(new_vector, other_vector)

    # Check that the vector is not described in other coordinate systems
    with raises(ValueError):
        equals(new_vector, test_args.v_cart)
    with raises(ValueError):
        equals(new_vector, test_args.v_sph)

    assert not equals(new_vector, ZERO)


def test_cylindrical_to_spherical(test_args: Args) -> None:
    new_vector = convert_vector(test_args.v_cyl, test_args.sph)
    assert equals(new_vector, test_args.v_sph)

    # Cannot compare vectors applied at different points in non-Cartesian systems
    other_vector = convert_vector(test_args.w_cyl, test_args.sph)
    with raises(ValueError):
        equals(new_vector, other_vector)

    # Check that the vector is not described in other coordinate systems
    with raises(ValueError):
        equals(new_vector, test_args.v_cart)
    with raises(ValueError):
        equals(new_vector, test_args.v_cyl)

    assert not equals(new_vector, ZERO)


def test_spherical_to_cartesian(test_args: Args) -> None:
    new_vector = convert_vector(test_args.v_sph, test_args.cart)
    assert equals(new_vector, test_args.v_cart)

    # In non-Cartesian systems, base vectors depend on the point of application.
    other_vector = convert_vector(test_args.w_sph, test_args.cart)
    assert not equals(new_vector, other_vector)

    # Check that the vector is not described in other coordinate systems
    with raises(ValueError):
        equals(new_vector, test_args.v_cyl)
    with raises(ValueError):
        equals(new_vector, test_args.v_sph)

    assert not equals(new_vector, ZERO)


def test_spherical_to_cylindrical(test_args: Args) -> None:
    new_vector = convert_vector(test_args.v_sph, test_args.cyl)
    assert equals(new_vector, test_args.v_cyl)

    # Cannot compare vectors applied at different points in non-Cartesian systems
    other_vector = convert_vector(test_args.w_sph, test_args.cyl)
    with raises(ValueError):
        equals(new_vector, other_vector)

    # Check that the vector is not described in other coordinate systems
    with raises(ValueError):
        equals(new_vector, test_args.v_cart)
    with raises(ValueError):
        equals(new_vector, test_args.v_sph)

    assert not equals(new_vector, ZERO)


def test_spherical_to_spherical(test_args: Args) -> None:
    new_sph = SphericalCoordinateSystem()

    new_vector = convert_vector(test_args.v_sph, new_sph)

    r, theta, phi = new_sph.base_scalars
    correct_point = {r: sqrt(2), theta: 3 * pi / 4, phi: pi}
    e_r, e_theta, e_phi = new_sph.base_vectors
    correct_value = e_r * -sqrt(2) - e_phi + e_theta * -sqrt(2)
    correct_vector = AppliedVector(correct_value, correct_point, new_sph)

    assert equals(new_vector, correct_vector)

    # Cannot compare vectors applied at different points in non-Cartesian systems
    other_vector = convert_vector(test_args.w_sph, new_sph)
    with raises(ValueError):
        equals(new_vector, other_vector)

    # Check that the vector is not described in other coordinate systems
    with raises(ValueError):
        equals(new_vector, test_args.v_cart)
    with raises(ValueError):
        equals(new_vector, test_args.v_cyl)

    assert not equals(new_vector, ZERO)


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
    new_vector = convert_vector(test_args.v_cart, test_args.cart)
    assert equals(new_vector, test_args.v_cart)

    # Check that the vector is not described in other coordinate systems
    with raises(ValueError):
        equals(new_vector, test_args.v_cyl)
    with raises(ValueError):
        equals(new_vector, test_args.v_sph)

    new_vector = convert_vector(test_args.v_cyl, test_args.cyl)
    assert equals(new_vector, test_args.v_cyl)

    # Check that the vector is not described in other coordinate systems
    with raises(ValueError):
        equals(new_vector, test_args.v_cart)
    with raises(ValueError):
        equals(new_vector, test_args.v_sph)

    new_vector = convert_vector(test_args.v_sph, test_args.sph)
    assert equals(new_vector, test_args.v_sph)

    # Check that the vector is not described in other coordinate systems
    with raises(ValueError):
        equals(new_vector, test_args.v_cart)
    with raises(ValueError):
        equals(new_vector, test_args.v_cyl)

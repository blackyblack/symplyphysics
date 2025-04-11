from dataclasses import dataclass
from pytest import fixture, raises
from sympy import pi, sqrt
from symplyphysics.core.experimental.solvers import vector_equals
from symplyphysics.core.experimental.vectors import VectorExpr
from symplyphysics.core.experimental.points import AppliedPoint
from symplyphysics.core.experimental.coordinate_systems import (
    CartesianCoordinateSystem,
    CylindricalCoordinateSystem,
    SphericalCoordinateSystem,
    express_base_vectors,
    convert_vector,
)


@dataclass(frozen=True)
class Args:  # pylint: disable=too-many-instance-attributes
    cart: CartesianCoordinateSystem
    cyl: CylindricalCoordinateSystem
    sph: SphericalCoordinateSystem

    p_cart: AppliedPoint
    p_cyl: AppliedPoint
    p_sph: AppliedPoint

    # These represent the same vector `v` in different coordinate systems
    v_cart: VectorExpr
    v_cyl: VectorExpr
    v_sph: VectorExpr

    q_cart: AppliedPoint
    q_cyl: AppliedPoint
    q_sph: AppliedPoint

    # These vectors have the same value as the respective `v` vector, but are applied at at
    # different points; they're not necessarily all equal to each other
    w_cart: VectorExpr
    w_cyl: VectorExpr
    w_sph: VectorExpr


@fixture(name="test_args")
def test_args_fixture() -> Args:
    cart = CartesianCoordinateSystem()
    cyl = CylindricalCoordinateSystem()
    sph = SphericalCoordinateSystem()

    # These represent the same physical point `p` (which is different from `q`) in different
    # coordinate systems
    p_cart = AppliedPoint([-1, 0, -1], cart)
    p_cyl = AppliedPoint([1, pi, -1], cyl)
    p_sph = AppliedPoint([sqrt(2), 3 * pi / 4, pi], sph)

    v_cart = cart.j + cart.k * 2

    _, e_phi_cyl_p, e_z_p = cyl.base_vectors(p_cyl)
    v_cyl = -e_phi_cyl_p + e_z_p * 2

    e_r_p, e_theta_p, e_phi_sph_p = sph.base_vectors(p_sph)
    v_sph = e_r_p * -sqrt(2) - e_phi_sph_p + e_theta_p * -sqrt(2)

    # These represent the same physical point `q` (which is different from `p`) in different
    # coordinate systems
    q_cart = AppliedPoint([1, 1, 1], cart)
    q_cyl = AppliedPoint([sqrt(2), pi / 4, 1], cyl)
    q_sph = AppliedPoint([sqrt(3), pi / 4, pi / 4], sph)

    w_cart = cart.j + cart.k * 2

    _, e_phi_cyl_q, e_z_q = cyl.base_vectors(q_cyl)
    w_cyl = -e_phi_cyl_q + e_z_q * 2

    e_r_q, e_theta_q, e_phi_sph_q = sph.base_vectors(q_sph)
    w_sph = e_r_q * -sqrt(2) - e_phi_sph_q + e_theta_q * -sqrt(2)

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
        q_cart=q_cart,
        q_cyl=q_cyl,
        q_sph=q_sph,
        w_cart=w_cart,
        w_cyl=w_cyl,
        w_sph=w_sph,
    )


def test_cartesian_to_cartesian(test_args: Args) -> None:
    new_cart = CartesianCoordinateSystem()

    new_vector = convert_vector(test_args.v_cart, test_args.p_cart, new_cart)

    _, j, k = new_cart.base_vectors()
    correct_vector = j + k * 2

    assert vector_equals(new_vector, correct_vector)

    # The transformation of basis vectors between Cartesian coordinates does *not* depend how we
    # choose the point where the vector is applied.
    other_vector = convert_vector(test_args.w_cart, test_args.q_cart, new_cart)
    assert vector_equals(new_vector, other_vector)

    # Check that the vector is not described in other coordinate systems
    assert not vector_equals(new_vector, test_args.v_cyl)
    assert not vector_equals(new_vector, test_args.v_sph)

    assert not vector_equals(new_vector, 0)


def test_cartesian_to_cylindrical(test_args: Args) -> None:
    new_vector = convert_vector(test_args.v_cart, test_args.p_cart, test_args.cyl)
    assert vector_equals(new_vector, test_args.v_cyl)

    # Cannot compare vectors applied at different points in non-Cartesian systems
    other_vector = convert_vector(test_args.w_cart, test_args.q_cart, test_args.cyl)
    assert not vector_equals(new_vector, other_vector)

    # Check that the vector is not described in other coordinate systems
    assert not vector_equals(new_vector, test_args.v_cart)
    assert not vector_equals(new_vector, test_args.v_sph)

    assert not vector_equals(new_vector, 0)


def test_cartesian_to_spherical(test_args: Args) -> None:
    new_vector = convert_vector(test_args.v_cart, test_args.p_cart, test_args.sph)
    assert vector_equals(new_vector, test_args.v_sph)

    # Cannot compare vectors applied at different points in non-Cartesian systems
    other_vector = convert_vector(test_args.w_cart, test_args.q_cart, test_args.sph)
    assert not vector_equals(new_vector, other_vector)

    # Check that the vector is not described in other coordinate systems
    assert not vector_equals(new_vector, test_args.v_cart)
    assert not vector_equals(new_vector, test_args.v_cyl)

    assert not vector_equals(new_vector, 0)


def test_cylindrical_to_cartesian(test_args: Args) -> None:
    new_vector = convert_vector(test_args.v_cyl, test_args.p_cyl, test_args.cart)
    assert vector_equals(new_vector, test_args.v_cart)

    # In non-Cartesian systems, base vectors depend on the point of application.
    other_vector = convert_vector(test_args.w_cyl, test_args.q_cyl, test_args.cart)
    assert not vector_equals(new_vector, other_vector)

    # Check that the vector is not described in other coordinate systems
    assert not vector_equals(new_vector, test_args.v_cyl)
    assert not vector_equals(new_vector, test_args.v_sph)

    assert not vector_equals(new_vector, 0)


def test_cylindrical_to_cylindrical(test_args: Args) -> None:
    new_cyl = CylindricalCoordinateSystem()

    new_vector = convert_vector(test_args.v_cyl, test_args.p_cyl, new_cyl)

    correct_point = AppliedPoint([1, pi, -1], new_cyl)
    _, e_phi, e_z = new_cyl.base_vectors(correct_point)
    correct_vector = -e_phi + e_z * 2

    assert vector_equals(new_vector, correct_vector)

    # Cannot compare vectors applied at different points in non-Cartesian systems
    other_vector = convert_vector(test_args.w_cyl, test_args.q_cyl, new_cyl)
    assert not vector_equals(new_vector, other_vector)

    # Check that the vector is not described in other coordinate systems
    assert not vector_equals(new_vector, test_args.v_cart)
    assert not vector_equals(new_vector, test_args.v_sph)

    assert not vector_equals(new_vector, 0)


def test_cylindrical_to_spherical(test_args: Args) -> None:
    new_vector = convert_vector(test_args.v_cyl, test_args.p_cyl, test_args.sph)
    assert vector_equals(new_vector, test_args.v_sph)

    # Cannot compare vectors applied at different points in non-Cartesian systems
    other_vector = convert_vector(test_args.w_cyl, test_args.q_cyl, test_args.sph)
    assert not vector_equals(new_vector, other_vector)

    # Check that the vector is not described in other coordinate systems
    assert not vector_equals(new_vector, test_args.v_cart)
    assert not vector_equals(new_vector, test_args.v_cyl)

    assert not vector_equals(new_vector, 0)


def test_spherical_to_cartesian(test_args: Args) -> None:
    new_vector = convert_vector(test_args.v_sph, test_args.p_sph, test_args.cart)
    assert vector_equals(new_vector, test_args.v_cart)

    # In non-Cartesian systems, base vectors depend on the point of application.
    other_vector = convert_vector(test_args.w_sph, test_args.q_sph, test_args.cart)
    assert not vector_equals(new_vector, other_vector)

    # Check that the vector is not described in other coordinate systems
    assert not vector_equals(new_vector, test_args.v_cyl)
    assert not vector_equals(new_vector, test_args.v_sph)

    assert not vector_equals(new_vector, 0)


def test_spherical_to_cylindrical(test_args: Args) -> None:
    new_vector = convert_vector(test_args.v_sph, test_args.p_sph, test_args.cyl)
    assert vector_equals(new_vector, test_args.v_cyl)

    # Cannot compare vectors applied at different points in non-Cartesian systems
    other_vector = convert_vector(test_args.w_sph, test_args.q_sph, test_args.cyl)
    assert not vector_equals(new_vector, other_vector)

    # Check that the vector is not described in other coordinate systems
    assert not vector_equals(new_vector, test_args.v_cart)
    assert not vector_equals(new_vector, test_args.v_sph)

    assert not vector_equals(new_vector, 0)


def test_spherical_to_spherical(test_args: Args) -> None:
    new_sph = SphericalCoordinateSystem()

    new_vector = convert_vector(test_args.v_sph, test_args.p_sph, new_sph)

    new_point = AppliedPoint([sqrt(2), 3 * pi / 4, pi], new_sph)
    e_r, e_theta, e_phi = new_sph.base_vectors(new_point)
    correct_vector = e_r * -sqrt(2) - e_phi + e_theta * -sqrt(2)

    assert vector_equals(new_vector, correct_vector)

    # Cannot compare vectors applied at different points in non-Cartesian systems
    other_vector = convert_vector(test_args.w_sph, test_args.q_sph, new_sph)
    assert not vector_equals(new_vector, other_vector)

    # Check that the vector is not described in other coordinate systems
    assert not vector_equals(new_vector, test_args.v_cart)
    assert not vector_equals(new_vector, test_args.v_cyl)

    assert not vector_equals(new_vector, 0)


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
    new_vector = convert_vector(test_args.v_cart, test_args.p_cart, test_args.cart)
    assert vector_equals(new_vector, test_args.v_cart)

    # Check that the vector is not described in other coordinate systems
    assert not vector_equals(new_vector, test_args.v_cyl)
    assert not vector_equals(new_vector, test_args.v_sph)

    new_vector = convert_vector(test_args.v_cyl, test_args.p_cyl, test_args.cyl)
    assert vector_equals(new_vector, test_args.v_cyl)

    # Check that the vector is not described in other coordinate systems
    assert not vector_equals(new_vector, test_args.v_cart)
    assert not vector_equals(new_vector, test_args.v_sph)

    new_vector = convert_vector(test_args.v_sph, test_args.p_sph, test_args.sph)
    assert vector_equals(new_vector, test_args.v_sph)

    # Check that the vector is not described in other coordinate systems
    assert not vector_equals(new_vector, test_args.v_cart)
    assert not vector_equals(new_vector, test_args.v_cyl)

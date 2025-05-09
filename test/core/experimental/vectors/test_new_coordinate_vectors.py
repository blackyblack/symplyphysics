from dataclasses import dataclass
from pytest import raises, fixture
from sympy import Symbol as SymSymbol, Function as SymFunction, ImmutableMatrix, sin, cos, sqrt, symbols as sym_symbols
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.experimental.vectors import (is_vector_expr, VectorSymbol, VectorNorm as
    norm, VectorDot as dot, VectorCross as cross, vector_diff)
from symplyphysics.core.experimental.coordinate_systems.new_coordinate_systems import (
    CartesianCoordinateSystem, CylindricalCoordinateSystem, SphericalCoordinateSystem)
from symplyphysics.core.experimental.vectors.new_coordinate_vectors import CoordinateVector, combine_coordinate_vectors


@dataclass(frozen=True, kw_only=True)
class Args:
    cart: CartesianCoordinateSystem
    cyl: CylindricalCoordinateSystem
    sph: SphericalCoordinateSystem
    t: SymSymbol


@fixture(name="test_args")
def test_args_fixture() -> Args:
    cart = CartesianCoordinateSystem()
    cyl = CylindricalCoordinateSystem()
    sph = SphericalCoordinateSystem()
    t = SymSymbol("t", real=True)

    return Args(cart=cart, cyl=cyl, sph=sph, t=t)


def test_cartesian_coordinate_vector(test_args: Args) -> None:
    v1 = CoordinateVector([1, 1, 0], test_args.cart)
    assert v1.components == ImmutableMatrix([1, 1, 0])
    assert v1.system is test_args.cart

    v2 = CoordinateVector([-1, -1, -1], test_args.cart)
    assert v2.system is test_args.cart
    assert v1.system is v2.system

    assert is_vector_expr(v1)
    assert is_vector_expr(v2 * 2)
    assert is_vector_expr(v1 + v2)

    assert combine_coordinate_vectors(v1 + v2) == CoordinateVector([0, 0, -1], test_args.cart)

    assert combine_coordinate_vectors(v1 * 0) == 0

    assert dot(v1, v2) == -2
    assert norm(v1) == sqrt(2)

    v3 = cross(v1, v2)
    assert isinstance(v3, CoordinateVector)
    assert v3 == CoordinateVector([-1, 1, 0], test_args.cart)

    t = test_args.t
    assert vector_diff(v1, t) == 0

    vt = CoordinateVector([t, t, -t], test_args.cart)
    assert vector_diff(vt) == CoordinateVector([1, 1, -1], test_args.cart)


def test_cylindrical_coordinate_vector(test_args: Args) -> None:
    v1 = CoordinateVector([1, 1, 1], test_args.cyl)
    assert v1.components == ImmutableMatrix([1, 1, 1])
    assert v1.system is test_args.cyl

    v2 = CoordinateVector([2, -1, -2], test_args.cyl)
    assert v2.system is v1.system

    assert is_vector_expr(v1)
    assert is_vector_expr(v2)
    assert is_vector_expr(v1 + v2)
    assert is_vector_expr(-v1 + 3 * v2)

    assert combine_coordinate_vectors(v1 - v2) == CoordinateVector([-1, 2, 3], test_args.cyl)

    assert dot(v1, v2) == -1
    assert norm(v2) == 3

    v3 = cross(v1, v2)
    assert v3 == CoordinateVector([-1, 4, -3], test_args.cyl)

    _, phi_f, _ = test_args.cyl.base_scalar_functions
    g = sym_symbols("g", cls=SymFunction)
    t = test_args.t
    v = CoordinateVector([g(t) + 1, 1, 1], test_args.cyl)
    dv_dt = CoordinateVector([-2 * t, 2 * t, 0], test_args.cyl)
    assert vector_diff(v, t).subs(g(t), 0).replace(phi_f, lambda t: 1 + t**2).doit() == dv_dt


def test_spherical_coordinate_vector(test_args: Args) -> None:
    v1 = CoordinateVector([1, 1, 1], test_args.sph)
    assert v1.components == ImmutableMatrix([1, 1, 1])
    assert v1.system is test_args.sph

    v2 = CoordinateVector([2, -1, -2], test_args.sph)
    assert v2.system is v1.system

    assert is_vector_expr(v1)
    assert is_vector_expr(v2)
    assert is_vector_expr(v1 + v2)
    assert is_vector_expr(-v1 + 3 * v2)

    assert combine_coordinate_vectors(v1 - v2) == CoordinateVector([-1, 2, 3], test_args.sph)

    assert dot(v1, v2) == -1
    assert norm(combine_coordinate_vectors(v1 - v2)) == sqrt(14)

    v3 = cross(v1, v2)
    assert v3 == CoordinateVector([-1, 4, -3], test_args.sph)

    _, theta_f, phi_f = test_args.sph.base_scalar_functions
    g = sym_symbols("g", cls=SymFunction)
    t = test_args.t
    v = CoordinateVector([g(t) + 1, 1, 1], test_args.sph)
    dv_dt = CoordinateVector(
        [-1 + sin(t) / t**2, 1 + cos(t) / t**2, -(sin(t) + cos(t)) / t**2],
        test_args.sph,
    )
    actual = vector_diff(v, t).subs(g(t), 0).replace(theta_f,
        lambda t: t).replace(phi_f, lambda t: 1 / t).doit()
    assert isinstance(actual, CoordinateVector)

    for lhs, rhs in zip(actual.components, dv_dt.components):
        assert expr_equals(lhs, rhs)

    assert actual.system is dv_dt.system


def test_wrong_components(test_args: Args) -> None:
    components = ImmutableMatrix([[1, 1], [2, 2]])

    with raises(ValueError):
        _ = CoordinateVector(components, test_args.cyl)

    components = [1, 2, "3"]
    with raises(ValueError):
        _ = CoordinateVector(components, test_args.sph)


def test_combine_coordinate_vectors(test_args: Args) -> None:
    v1 = CoordinateVector([0, 0, 0], test_args.cart)
    v2 = CoordinateVector([1, -1, 1], test_args.cart)
    v3 = CoordinateVector([1, 1, 1], test_args.cart)

    v4 = CoordinateVector([-1, -1, -1], test_args.cyl)
    v5 = CoordinateVector([1, 1, 1], test_args.cyl)

    v6 = CoordinateVector([4, 4, 4], test_args.sph)

    v_coord = v1 + v2 + v3 + v4 + v5 + v6

    v7 = VectorSymbol("v_7")
    v8 = VectorSymbol("v_8")
    v9 = cross(v7, v8)

    v_symbo = v7 + v8 + v9

    v = v_coord + v_symbo

    actual = combine_coordinate_vectors(v)

    combined = (CoordinateVector([2, 0, 2], test_args.cart) +
        CoordinateVector([4, 4, 4], test_args.sph) + v7 + v8 + v9)

    assert norm(combine_coordinate_vectors(actual - combined)) == 0

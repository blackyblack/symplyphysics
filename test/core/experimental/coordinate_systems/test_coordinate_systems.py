from pytest import raises
from sympy import sin
from symplyphysics import Symbol, units, angle_type
from symplyphysics.core.errors import UnitsError
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.experimental.vectors import (VectorSymbol, VectorFunction, VectorNorm as
    norm, VectorDot as dot, VectorCross as cross, AppliedVectorFunction)
from symplyphysics.core.experimental.coordinate_systems import (
    CartesianCoordinateSystem,
    CylindricalCoordinateSystem,
    SphericalCoordinateSystem,
)
from symplyphysics.core.experimental.points import AppliedPoint


def test_base_init() -> None:
    x = Symbol("x", units.length)
    y = Symbol("y", units.length)
    z = Symbol("z", units.length)
    i = VectorSymbol("i")
    j = VectorSymbol("j")
    k = VectorSymbol("k")

    good_scalars = [x, y, z]
    good_vectors = [i, j, k]

    # If the base scalars are specified, there should be exactly 3 of them
    for idx in range(1, 3):
        with raises(ValueError):
            CartesianCoordinateSystem(base_scalars=good_scalars[:idx])

    # If the base vectors are specified, there should be exactly 3 of them
    for idx in range(1, 3):
        with raises(ValueError):
            CartesianCoordinateSystem(base_vectors=good_vectors[:idx])

    # The base vectors must be unique
    with raises(ValueError):
        CartesianCoordinateSystem(base_vectors=[i, i, k])
    with raises(ValueError):
        CartesianCoordinateSystem(base_vectors=[j, j, j])


def test_cartesian_system() -> None:
    # __init__ without arguments

    c = CartesianCoordinateSystem()

    assert len(c.base_scalars) == 3
    assert c.x == c.base_scalars[0]
    assert c.y == c.base_scalars[1]
    assert c.z == c.base_scalars[2]

    assert len(c.base_vectors()) == 3
    assert c.i, c.base_vectors()[0]
    assert c.j == c.base_vectors()[1]
    assert c.k == c.base_vectors()[2]

    assert not isinstance(c.i, AppliedVectorFunction)
    assert not isinstance(c.j, AppliedVectorFunction)
    assert not isinstance(c.k, AppliedVectorFunction)

    assert len(c.lame_coefficients) == 3
    assert all(h == 1 for h in c.lame_coefficients)

    assert c.jacobian == 1

    # __init__ with arguments

    x = Symbol("x", units.length)
    y = Symbol("y", units.length)
    z = Symbol("z", units.length)
    i = VectorSymbol("i")
    j = VectorSymbol("j")
    k = VectorSymbol("k")

    c = CartesianCoordinateSystem(
        base_scalars=[x, y, z],
        base_vectors=[i, j, k],
    )

    assert c.x == x
    assert c.y == y
    assert c.z == z
    assert c.i == i
    assert c.j == j
    assert c.k == k

    # errors

    # Wrong dimension of base scalars
    good_scalars = [x, y, z]
    bad_scalar = Symbol("s", units.area)
    for idx in range(3):
        bad_scalars = good_scalars.copy()
        bad_scalars[idx] = bad_scalar

        # NOTE: for now, making `bad_scalar` dimensionless would result in `TypeError`, not `UnitsError`
        with raises(UnitsError):
            CartesianCoordinateSystem(base_scalars=bad_scalars)

    # Wrong dimension of base vectors
    good_vectors = [i, j, k]
    bad_vector = VectorSymbol("i", units.length)
    for idx in range(3):
        bad_vectors = good_vectors.copy()
        bad_vectors[idx] = bad_vector
        with raises(UnitsError):
            CartesianCoordinateSystem(base_vectors=bad_vectors)

    # expr_simplify

    kvs = {
        norm(i): 1,
        norm(j): 1,
        norm(k): 1,
        dot(i, i): 1,
        dot(i, j): 0,
        dot(i, k): 0,
        dot(j, i): 0,
        dot(j, j): 1,
        dot(j, k): 0,
        dot(k, i): 0,
        dot(k, j): 0,
        dot(k, k): 1,
        cross(i, i): 0,
        cross(i, j): k,
        cross(i, k): -j,
        cross(j, i): -k,
        cross(j, j): 0,
        cross(j, k): i,
        cross(k, i): j,
        cross(k, j): -i,
        cross(k, k): 0,
    }

    for expr, expected in kvs.items():
        assert expr_equals(c.expr_simplify(expr), expected)


def test_cylindrical_system() -> None:
    k = CartesianCoordinateSystem()
    p = AppliedPoint(k.base_scalars, k)

    # __init__ without arguments

    c = CylindricalCoordinateSystem()

    assert len(c.base_scalars) == 3
    assert c.rho == c.base_scalars[0]
    assert c.phi == c.base_scalars[1]
    assert c.z == c.base_scalars[2]

    assert len(c.base_vectors(p)) == 3
    e_rho, e_phi, e_z = c.base_vectors(p)

    assert isinstance(e_rho, AppliedVectorFunction)
    assert isinstance(e_phi, AppliedVectorFunction)
    assert not isinstance(e_z, AppliedVectorFunction)

    assert len(c.lame_coefficients) == 3
    assert c.lame_coefficients[0] == 1
    assert c.lame_coefficients[1] == c.rho
    assert c.lame_coefficients[2] == 1

    assert expr_equals(c.jacobian, c.rho)

    # __init__ with arguments

    rho = Symbol("rho", units.length)
    phi = Symbol("phi")
    z = Symbol("z", units.length)
    e_rho = VectorFunction("e_rho", arguments=(p,))
    e_phi = VectorFunction("e_phi", arguments=(p,))
    e_z = VectorSymbol("e_z")

    c = CylindricalCoordinateSystem(
        base_scalars=[rho, phi, z],
        base_vectors=[e_rho, e_phi, e_z],
    )

    assert c.rho == rho
    assert c.phi == phi
    assert c.z == z
    assert c.base_vectors(p)[0] == e_rho(p)
    assert c.base_vectors(p)[1] == e_phi(p)
    assert c.base_vectors(p)[2] == e_z

    # errors

    # Wrong dimension of base scalars
    good_scalars = [rho, phi, z]
    bad_scalar = Symbol("S", units.energy / units.temperature)
    for idx in range(3):
        bad_scalars = good_scalars.copy()
        bad_scalars[idx] = bad_scalar

        # NOTE: for now, making `bad_scalar` dimensionless would result in `TypeError`, not `UnitsError`
        with raises(UnitsError):
            CylindricalCoordinateSystem(base_scalars=bad_scalars)

    # Wrong dimension of base vectors
    good_vectors: list[VectorSymbol | VectorFunction] = [e_rho, e_phi, e_z]
    bad_vector = VectorSymbol("a", units.force)
    for idx in range(3):
        bad_vectors = good_vectors.copy()
        bad_vectors[idx] = bad_vector
        with raises(UnitsError):
            CylindricalCoordinateSystem(base_vectors=bad_vectors)

    # expr_simplify

    kvs = {
        norm(e_rho(p)): 1,
        norm(e_phi(p)): 1,
        norm(e_z): 1,
        dot(e_rho(p), e_rho(p)): 1,
        dot(e_rho(p), e_phi(p)): 0,
        dot(e_rho(p), e_z): 0,
        dot(e_phi(p), e_rho(p)): 0,
        dot(e_phi(p), e_phi(p)): 1,
        dot(e_phi(p), e_z): 0,
        dot(e_z, e_rho(p)): 0,
        dot(e_z, e_phi(p)): 0,
        dot(e_z, e_z): 1,
        cross(e_rho(p), e_rho(p)): 0,
        cross(e_rho(p), e_phi(p)): e_z,
        cross(e_rho(p), e_z): -e_phi(p),
        cross(e_phi(p), e_rho(p)): -e_z,
        cross(e_phi(p), e_phi(p)): 0,
        cross(e_phi(p), e_z): e_rho(p),
        cross(e_z, e_rho(p)): e_phi(p),
        cross(e_z, e_phi(p)): -e_rho(p),
        cross(e_z, e_z): 0,
    }

    for expr, expected in kvs.items():
        assert expr_equals(c.expr_simplify(expr, p), expected)


def test_spherical_system() -> None:
    k = CartesianCoordinateSystem()
    p = AppliedPoint(k.base_scalars, k)

    # __init__ without arguments

    c = SphericalCoordinateSystem()

    assert len(c.base_scalars) == 3
    assert c.r == c.base_scalars[0]
    assert c.theta == c.base_scalars[1]
    assert c.phi == c.base_scalars[2]

    assert len(c.base_vectors(p)) == 3
    e_r, e_theta, e_phi = c.base_vectors(p)
    assert isinstance(e_r, AppliedVectorFunction)
    assert isinstance(e_theta, AppliedVectorFunction)
    assert isinstance(e_phi, AppliedVectorFunction)

    assert len(c.lame_coefficients) == 3
    assert c.lame_coefficients[0] == 1
    assert c.lame_coefficients[1] == c.r
    assert c.lame_coefficients[2] == c.r * sin(c.theta)

    assert expr_equals(c.jacobian, c.r**2 * sin(c.theta))

    # __init__ with arguments

    r = Symbol("r", units.length)
    theta = Symbol("theta", angle_type)
    phi = Symbol("phi")
    e_r = VectorFunction("e_r", arguments=(p,))
    e_theta = VectorFunction("e_theta", arguments=(p,))
    e_phi = VectorFunction("e_phi", arguments=(p,))

    c = SphericalCoordinateSystem(
        base_scalars=[r, theta, phi],
        base_vectors=[e_r, e_theta, e_phi],
    )

    assert c.r == r
    assert c.theta == theta
    assert c.phi == phi
    assert c.base_vectors(p)[0] == e_r(p)
    assert c.base_vectors(p)[1] == e_theta(p)
    assert c.base_vectors(p)[2] == e_phi(p)

    # errors

    # Wrong dimension of base scalars
    good_scalars = [r, theta, phi]
    bad_scalar = Symbol("T", units.temperature)
    for idx in range(3):
        bad_scalars = good_scalars.copy()
        bad_scalars[idx] = bad_scalar

        # NOTE: for now, making `bad_scalar` dimensionless would result in `TypeError`, not `UnitsError`
        with raises(UnitsError):
            SphericalCoordinateSystem(base_scalars=bad_scalars)

    # Wrong dimension of base vectors
    good_vectors = [e_r, e_theta, e_phi]
    bad_vector = VectorSymbol("E", units.force / units.charge)
    for idx in range(3):
        bad_vectors = good_vectors.copy()
        bad_vectors[idx] = bad_vector
        with raises(UnitsError):
            SphericalCoordinateSystem(base_vectors=bad_vectors)

    # expr_simplify

    kvs = {
        norm(e_r(p)): 1,
        norm(e_theta(p)): 1,
        norm(e_phi(p)): 1,
        dot(e_r(p), e_r(p)): 1,
        dot(e_r(p), e_theta(p)): 0,
        dot(e_r(p), e_phi(p)): 0,
        dot(e_theta(p), e_r(p)): 0,
        dot(e_theta(p), e_theta(p)): 1,
        dot(e_theta(p), e_phi(p)): 0,
        dot(e_phi(p), e_r(p)): 0,
        dot(e_phi(p), e_theta(p)): 0,
        dot(e_phi(p), e_phi(p)): 1,
        cross(e_r(p), e_r(p)): 0,
        cross(e_r(p), e_theta(p)): e_phi(p),
        cross(e_r(p), e_phi(p)): -e_theta(p),
        cross(e_theta(p), e_r(p)): -e_phi(p),
        cross(e_theta(p), e_theta(p)): 0,
        cross(e_theta(p), e_phi(p)): e_r(p),
        cross(e_phi(p), e_r(p)): e_theta(p),
        cross(e_phi(p), e_theta(p)): -e_r(p),
        cross(e_phi(p), e_phi(p)): 0,
    }

    for expr, expected in kvs.items():
        assert expr_equals(c.expr_simplify(expr, p), expected)

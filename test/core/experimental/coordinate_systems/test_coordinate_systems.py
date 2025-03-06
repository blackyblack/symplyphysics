from pytest import raises
from sympy import sin
from symplyphysics import Symbol, units, angle_type
from symplyphysics.core.errors import UnitsError
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.experimental.vectors import VectorSymbol
from symplyphysics.core.experimental.coordinate_systems import (
    CartesianCoordinateSystem,
    CylindricalCoordinateSystem,
    SphericalCoordinateSystem,
)


def test_base_init() -> None:
    x = Symbol("x", units.length)
    y = Symbol("y", units.length)
    z = Symbol("z", units.length)
    i = VectorSymbol("i", norm=1)
    j = VectorSymbol("j", norm=1)
    k = VectorSymbol("k", norm=1)

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

    assert len(c.base_vectors) == 3
    assert c.i == c.base_vectors[0]
    assert c.j == c.base_vectors[1]
    assert c.k == c.base_vectors[2]

    assert len(c.lame_coefficients) == 3
    assert all(h == 1 for h in c.lame_coefficients)

    assert c.jacobian == 1

    # __init__ with arguments

    x = Symbol("x", units.length)
    y = Symbol("y", units.length)
    z = Symbol("z", units.length)
    i = VectorSymbol("i", norm=1)
    j = VectorSymbol("j", norm=1)
    k = VectorSymbol("k", norm=1)

    c = CartesianCoordinateSystem(base_scalars=[x, y, z], base_vectors=[i, j, k])

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


def test_cylindrical_system() -> None:
    # __init__ without arguments

    c = CylindricalCoordinateSystem()

    assert len(c.base_scalars) == 3
    assert c.rho == c.base_scalars[0]
    assert c.phi == c.base_scalars[1]
    assert c.z == c.base_scalars[2]

    assert len(c.base_vectors) == 3

    assert len(c.lame_coefficients) == 3
    assert c.lame_coefficients[0] == 1
    assert c.lame_coefficients[1] == c.rho
    assert c.lame_coefficients[2] == 1

    assert expr_equals(c.jacobian, c.rho)

    # __init__ with arguments

    rho = Symbol("rho", units.length)
    phi = Symbol("phi")
    z = Symbol("z", units.length)
    e_rho = VectorSymbol("e_rho", norm=1)
    e_phi = VectorSymbol("e_phi", norm=1)
    e_z = VectorSymbol("e_z", norm=1)

    c = CylindricalCoordinateSystem(base_scalars=[rho, phi, z], base_vectors=[e_rho, e_phi, e_z])

    assert c.rho == rho
    assert c.phi == phi
    assert c.z == z
    assert c.base_vectors[0] == e_rho
    assert c.base_vectors[1] == e_phi
    assert c.base_vectors[2] == e_z

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
    good_vectors = [e_rho, e_phi, e_z]
    bad_vector = VectorSymbol("a", units.force)
    for idx in range(3):
        bad_vectors = good_vectors.copy()
        bad_vectors[idx] = bad_vector
        with raises(UnitsError):
            CylindricalCoordinateSystem(base_vectors=bad_vectors)


def test_spherical_system() -> None:
    # __init__ without arguments

    c = SphericalCoordinateSystem()

    assert len(c.base_scalars) == 3
    assert c.r == c.base_scalars[0]
    assert c.theta == c.base_scalars[1]
    assert c.phi == c.base_scalars[2]

    assert len(c.base_vectors) == 3

    assert len(c.lame_coefficients) == 3
    assert c.lame_coefficients[0] == 1
    assert c.lame_coefficients[1] == c.r
    assert c.lame_coefficients[2] == c.r * sin(c.theta)

    assert expr_equals(c.jacobian, c.r**2 * sin(c.theta))

    # __init__ with arguments

    r = Symbol("r", units.length)
    theta = Symbol("theta", angle_type)
    phi = Symbol("phi")
    e_rho = VectorSymbol("e_rho", norm=1)
    e_theta = VectorSymbol("e_theta", norm=1)
    e_phi = VectorSymbol("e_phi", norm=1)

    c = SphericalCoordinateSystem(
        base_scalars=[r, theta, phi],
        base_vectors=[e_rho, e_theta, e_phi],
    )

    assert c.r == r
    assert c.theta == theta
    assert c.phi == phi
    assert c.base_vectors[0] == e_rho
    assert c.base_vectors[1] == e_theta
    assert c.base_vectors[2] == e_phi

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
    good_vectors = [e_rho, e_theta, e_phi]
    bad_vector = VectorSymbol("E", units.force / units.charge)
    for idx in range(3):
        bad_vectors = good_vectors.copy()
        bad_vectors[idx] = bad_vector
        with raises(UnitsError):
            SphericalCoordinateSystem(base_vectors=bad_vectors)

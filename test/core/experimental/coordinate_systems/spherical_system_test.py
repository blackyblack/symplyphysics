from dataclasses import dataclass
from pytest import fixture
from sympy import ImmutableMatrix, sin, cos, atan2, sqrt, pi, S
from symplyphysics import Symbol, units, clone_as_function
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.experimental.coordinate_systems import (BaseCoordinateSystem,
    SphericalCoordinateSystem)


@dataclass(kw_only=True, frozen=True)
# pylint: disable=too-many-instance-attributes
class Args:
    t: Symbol

    r: Symbol
    phi: Symbol
    theta: Symbol


@fixture(name="test_args")
def test_args_fixture() -> Args:
    t = Symbol("t", units.time)

    r = Symbol("r", units.length, positive=True)
    phi = Symbol("phi", real=True)
    theta = Symbol("theta", real=True)

    return Args(
        t=t,
        r=r,
        phi=phi,
        theta=theta,
    )


def test_spherical_matrix(test_args: Args) -> None:
    sph = SphericalCoordinateSystem((test_args.r, test_args.theta, test_args.phi))
    assert sph.base_scalars == (test_args.r, test_args.theta, test_args.phi)

    sph = SphericalCoordinateSystem()
    r, theta, phi = sph.base_scalars
    assert isinstance(sph, BaseCoordinateSystem)

    assert sph.cartesian_transform([1, pi / 3, pi / 2]) == (0, sqrt(3) / 2, S.Half)
    assert sph.cartesian_transform([2, pi, 0]) == (0, 0, -2)

    assert sph.inverse_cartesian_transform([1, 1, 1]) == (sqrt(3), atan2(sqrt(2), 1), pi / 4)
    assert sph.inverse_cartesian_transform([3, 4, 5]) == (5 * sqrt(2), pi / 4, atan2(4, 3))

    matrix = ImmutableMatrix([
        [0, sqrt(3) / 2, S.Half],
        [0, S.Half, -sqrt(3) / 2],
        [-sqrt(3) / 2, 0, 0],
    ])
    assert sph.cartesian_derivative_matrix([1, pi / 3, pi / 2]) == matrix

    assert sph.lame_coefficients([1, 5 * pi / 6, 3 * pi / 4]) == (1, 1, S.Half)

    matrix = ImmutableMatrix([
        [-sqrt(6) / 4, sqrt(6) / 4, S.Half],
        [-sqrt(2) / 4, sqrt(2) / 4, -sqrt(3) / 2],
        [-sqrt(2) / 2, -sqrt(2) / 2, 0],
    ])
    assert sph.base_vector_matrix([1, pi / 3, 3 * pi / 4]) == matrix
    assert matrix.det() == 1

    # Check that the system gives an orthogonal basis indeed
    mtx = sph.base_vector_matrix()
    v1, v2, v3 = (mtx.col(j) for j in range(3))
    assert expr_equals(v1.dot(v2), 0)
    assert expr_equals(v1.dot(v3), 0)
    assert expr_equals(v2.dot(v3), 0)

    t = test_args.t
    r_f = clone_as_function(r)
    theta_f = clone_as_function(theta)
    phi_f = clone_as_function(phi)

    dtheta_dt = theta_f(t).diff(t)
    dphi_dt = phi_f(t).diff(t)

    matrix = ImmutableMatrix([
        [0, dtheta_dt, sin(theta_f(t)) * dphi_dt],
        [-dtheta_dt, 0, cos(theta_f(t)) * dphi_dt],
        [-sin(theta_f(t)) * dphi_dt, -cos(theta_f(t)) * dphi_dt, 0],
    ])
    assert sph.diff_base_vector_matrix(t, [r_f(t), theta_f(t), phi_f(t)]) == matrix

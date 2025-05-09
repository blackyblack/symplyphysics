from typing import Sequence, Optional
from dataclasses import dataclass
from pytest import fixture, raises
from sympy import ImmutableMatrix, sin, cos, atan2, sqrt, pi, S, Expr
from symplyphysics import Symbol, units, clone_as_function
from symplyphysics.core.experimental.coordinate_systems.new_coordinate_systems import (
    BaseCoordinateSystem, CartesianCoordinateSystem, CylindricalCoordinateSystem,
    SphericalCoordinateSystem)


@dataclass(kw_only=True, frozen=True)
# pylint: disable=too-many-instance-attributes
class Args:
    t: Symbol

    x: Symbol
    y: Symbol
    z: Symbol

    rho: Symbol
    r: Symbol
    phi: Symbol
    theta: Symbol


@fixture(name="test_args")
def test_args_fixture() -> Args:
    t = Symbol("t", units.time)

    x = Symbol("x", units.length, real=True)
    y = Symbol("y", units.length, real=True)
    z = Symbol("z", units.length, real=True)

    rho = Symbol("rho", units.length, positive=True)
    r = Symbol("r", units.length, positive=True)
    phi = Symbol("phi", real=True)
    theta = Symbol("theta", real=True)

    return Args(t=t, x=x, y=y, z=z, rho=rho, r=r, phi=phi, theta=theta)


def test_cartesian_system(test_args: Args) -> None:
    cart = CartesianCoordinateSystem((test_args.x, test_args.y, test_args.z))
    assert cart.base_scalars == (test_args.x, test_args.y, test_args.z)

    cart = CartesianCoordinateSystem()
    assert len(cart.base_scalars) == 3
    assert isinstance(cart, BaseCoordinateSystem)

    assert cart.cartesian_transform() == cart.base_scalars
    assert cart.inverse_cartesian_transform(cart.base_scalars) == cart.base_scalars

    assert cart.cartesian_derivative_matrix() == ImmutableMatrix.eye(3)

    assert cart.lame_coefficients() == (1, 1, 1)

    assert cart.base_vector_matrix() == ImmutableMatrix.eye(3)

    assert cart.diff_base_vector_matrix(test_args.t) == ImmutableMatrix.zeros(3)


def test_cylindrical_system(test_args: Args) -> None:
    cyl = CylindricalCoordinateSystem((test_args.rho, test_args.phi, test_args.z))
    assert cyl.base_scalars == (test_args.rho, test_args.phi, test_args.z)

    cyl = CylindricalCoordinateSystem()
    rho, phi, z = cyl.base_scalars
    assert isinstance(cyl, BaseCoordinateSystem)

    assert cyl.cartesian_transform([2, pi / 2, 1]) == (0, 2, 1)
    assert cyl.cartesian_transform([sqrt(2), pi / 4, -1]) == (1, 1, -1)

    assert cyl.inverse_cartesian_transform([1, 1, 1]) == (sqrt(2), pi / 4, 1)
    assert cyl.inverse_cartesian_transform([3, 4, 5]) == (5, atan2(4, 3), 5)

    matrix = ImmutableMatrix([[S.Half, sqrt(3) / 2, 0], [-sqrt(3), 1, 0], [0, 0, 1]])
    assert cyl.cartesian_derivative_matrix([2, pi / 3, 1]) == matrix

    assert cyl.lame_coefficients([5, 1, 1]) == (1, 5, 1)

    matrix = ImmutableMatrix([[-1, 0, 0], [0, -1, 0], [0, 0, 1]])
    assert cyl.base_vector_matrix([3, pi, -1]) == matrix
    assert matrix.det() == 1

    t = test_args.t
    rho_f = clone_as_function(rho)
    phi_f = clone_as_function(phi)
    z_f = clone_as_function(z)
    matrix = ImmutableMatrix([[0, phi_f(t).diff(t), 0], [-phi_f(t).diff(t), 0, 0], [0, 0, 0]])

    assert cyl.diff_base_vector_matrix(t, [rho_f(t), phi_f(t), z_f(t)]) == matrix


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


def test_no_generate_base_scalars() -> None:

    # pylint: disable-next=abstract-method
    class NewCoordinateSystem(BaseCoordinateSystem):
        pass

    with raises(NotImplementedError):
        _ = NewCoordinateSystem()


def test_no_generate_lame_coefficients() -> None:

    # pylint: disable-next=abstract-method
    class NewCoordinateSystem(BaseCoordinateSystem):

        def generate_base_scalars(self) -> tuple[Symbol, Symbol, Symbol]:
            return (
                Symbol("a", units.length),
                Symbol("b", units.length),
                Symbol("c", units.length),
            )

        def cartesian_transform(
            self,
            base_scalars: Optional[Sequence[Expr]] = None,
        ) -> Sequence[Expr]:
            a, b, c = base_scalars or self.base_scalars

            return a + b, a - b, 2 * c

        def inverse_cartesian_transform(self, cartesian_scalars: Sequence[Expr]) -> Sequence[Expr]:
            x, y, z = cartesian_scalars

            return (x + y) / 2, (x - y) / 2, z / 2

    sys = NewCoordinateSystem()
    assert sys.lame_coefficients() == (sqrt(2), sqrt(2), 2)

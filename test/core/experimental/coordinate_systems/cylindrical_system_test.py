from dataclasses import dataclass
from pytest import fixture
from sympy import ImmutableMatrix, atan2, sqrt, pi, S
from symplyphysics import Symbol, units, clone_as_function
from symplyphysics.core.experimental.coordinate_systems import (BaseCoordinateSystem,
    CylindricalCoordinateSystem)


@dataclass(kw_only=True, frozen=True)
# pylint: disable=too-many-instance-attributes
class Args:
    t: Symbol

    rho: Symbol
    phi: Symbol
    z: Symbol


@fixture(name="test_args")
def test_args_fixture() -> Args:
    t = Symbol("t", units.time)

    rho = Symbol("rho", units.length, positive=True)
    phi = Symbol("phi", real=True)
    z = Symbol("z", units.length, real=True)

    return Args(t=t, z=z, rho=rho, phi=phi)


def test_cylindrical_system(test_args: Args) -> None:
    sys = CylindricalCoordinateSystem((test_args.rho, test_args.phi, test_args.z))
    assert sys.base_scalars == (test_args.rho, test_args.phi, test_args.z)

    sys = CylindricalCoordinateSystem()
    rho, phi, z = sys.base_scalars
    assert isinstance(sys, BaseCoordinateSystem)

    assert sys.cartesian_transform([2, pi / 2, 1]) == (0, 2, 1)
    assert sys.cartesian_transform([sqrt(2), pi / 4, -1]) == (1, 1, -1)

    assert sys.inverse_cartesian_transform([1, 1, 1]) == (sqrt(2), pi / 4, 1)
    assert sys.inverse_cartesian_transform([3, 4, 5]) == (5, atan2(4, 3), 5)

    matrix = ImmutableMatrix([[S.Half, sqrt(3) / 2, 0], [-sqrt(3), 1, 0], [0, 0, 1]])
    assert sys.cartesian_derivative_matrix([2, pi / 3, 1]) == matrix

    assert sys.lame_coefficients([5, 1, 1]) == (1, 5, 1)

    matrix = ImmutableMatrix([[-1, 0, 0], [0, -1, 0], [0, 0, 1]])
    assert sys.base_vector_matrix([3, pi, -1]) == matrix
    assert matrix.det() == 1

    # Check that the system gives an orthogonal basis indeed
    mtx = sys.base_vector_matrix()
    v1, v2, v3 = (mtx.col(j) for j in range(3))
    assert v1.dot(v2) == 0
    assert v1.dot(v3) == 0
    assert v2.dot(v3) == 0

    t = test_args.t
    rho_f = clone_as_function(rho)
    phi_f = clone_as_function(phi)
    z_f = clone_as_function(z)
    matrix = ImmutableMatrix([[0, phi_f(t).diff(t), 0], [-phi_f(t).diff(t), 0, 0], [0, 0, 0]])

    assert sys.diff_base_vector_matrix(t, [rho_f(t), phi_f(t), z_f(t)]) == matrix

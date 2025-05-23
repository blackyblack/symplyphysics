from typing import Sequence, Optional
from dataclasses import dataclass
from pytest import fixture
from sympy import ImmutableMatrix, sqrt, Expr
from symplyphysics import Symbol, units
from symplyphysics.core.experimental.coordinate_systems import (BaseCoordinateSystem,
    CartesianCoordinateSystem)


@dataclass(kw_only=True, frozen=True)
class Args:
    t: Symbol

    x: Symbol
    y: Symbol
    z: Symbol


@fixture(name="test_args")
def test_args_fixture() -> Args:
    t = Symbol("t", units.time)

    x = Symbol("x", units.length, real=True)
    y = Symbol("y", units.length, real=True)
    z = Symbol("z", units.length, real=True)

    return Args(t=t, x=x, y=y, z=z)


def test_cartesian_system(test_args: Args) -> None:
    sys = CartesianCoordinateSystem((test_args.x, test_args.y, test_args.z))
    assert sys.base_scalars == (test_args.x, test_args.y, test_args.z)

    sys = CartesianCoordinateSystem()
    assert len(sys.base_scalars) == 3
    assert isinstance(sys, BaseCoordinateSystem)

    assert sys.cartesian_transform() == sys.base_scalars
    assert sys.inverse_cartesian_transform(sys.base_scalars) == sys.base_scalars

    assert sys.cartesian_derivative_matrix() == ImmutableMatrix.eye(3)

    assert sys.lame_coefficients() == (1, 1, 1)

    # This also shows that the basis is orthogonal
    assert sys.base_vector_matrix() == ImmutableMatrix.eye(3)

    assert sys.diff_base_vector_matrix(test_args.t) == ImmutableMatrix.zeros(3)


def test_rotated_scaled_cartesian_system() -> None:

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

    assert sys.cartesian_transform([1, 2, -1]) == (3, -1, -2)
    assert sys.inverse_cartesian_transform([2, 2, -2]) == (2, 0, -1)

    matrix = ImmutableMatrix([[1, 1, 0], [1, -1, 0], [0, 0, 2]])
    assert sys.cartesian_derivative_matrix() == matrix

    assert sys.lame_coefficients() == (sqrt(2), sqrt(2), 2)

    matrix = ImmutableMatrix([
        [1 / sqrt(2), 1 / sqrt(2), 0],
        [1 / sqrt(2), -1 / sqrt(2), 0],
        [0, 0, 1],
    ])
    assert sys.base_vector_matrix() == matrix

    # Check that the system gives an orthogonal basis indeed
    mtx = sys.base_vector_matrix()
    v1, v2, v3 = (mtx.col(j) for j in range(3))
    assert v1.dot(v2) == 0
    assert v1.dot(v3) == 0
    assert v2.dot(v3) == 0

    t = Symbol("t")
    assert sys.diff_base_vector_matrix(t) == ImmutableMatrix.zeros(3)

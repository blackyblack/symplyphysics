from typing import Sequence, Optional
from pytest import raises
from sympy import ImmutableMatrix, sqrt, Expr, Rational
from symplyphysics import Symbol, units
from symplyphysics.core.experimental.coordinate_systems import BaseCoordinateSystem


def test_no_generate_base_scalars() -> None:

    # pylint: disable-next=abstract-method
    class NewCoordinateSystem(BaseCoordinateSystem):
        pass

    with raises(NotImplementedError):
        _ = NewCoordinateSystem()


def test_parabolic_system() -> None:

    class ParabolicCoordinateSystem(BaseCoordinateSystem):

        def generate_base_scalars(self) -> tuple[Symbol, Symbol, Symbol]:
            return (
                Symbol("sigma", sqrt(units.length), display_latex="\\sigma", real=True),
                Symbol("tau", sqrt(units.length), display_latex="\\tau", real=True),
                Symbol("z", units.length, real=True),
            )

        def cartesian_transform(
            self,
            base_scalars: Optional[Sequence[Expr]] = None,
        ) -> Sequence[Expr]:
            sigma, tau, z = base_scalars or self.base_scalars

            x = sigma * tau
            y = (tau**2 - sigma**2) / 2

            return x, y, z

        def inverse_cartesian_transform(self, cartesian_scalars: Sequence[Expr]) -> Sequence[Expr]:
            raise NotImplementedError("Several solutions exist")

    sys = ParabolicCoordinateSystem()

    assert sys.cartesian_transform([2, 4, 1]) == (8, 6, 1)

    with raises(NotImplementedError):
        _ = sys.inverse_cartesian_transform([1, 1, 1])

    assert sys.cartesian_derivative_matrix([3, 4, 5]) == ImmutableMatrix([
        [4, -3, 0],
        [3, 4, 0],
        [0, 0, 1],
    ])

    sigma, tau, _ = sys.base_scalars
    assert sys.lame_coefficients() == (sqrt(sigma**2 + tau**2), sqrt(sigma**2 + tau**2), 1)

    assert sys.base_vector_matrix([3, 4, 5]) == ImmutableMatrix([
        [Rational(4, 5), Rational(-3, 5), 0],
        [Rational(3, 5), Rational(4, 5), 0],
        [0, 0, 1],
    ])

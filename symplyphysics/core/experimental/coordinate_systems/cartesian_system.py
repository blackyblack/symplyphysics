from typing import Optional, Sequence

from sympy import Expr, Matrix, Symbol as SymSymbol, S
from sympy.physics import units

from symplyphysics.core.symbols.symbols import Symbol

from .base_system import BaseCoordinateSystem


class CartesianCoordinateSystem(BaseCoordinateSystem):

    def generate_base_scalars(self) -> tuple[Symbol, Symbol, Symbol]:
        return (
            Symbol("x", units.length, real=True),
            Symbol("y", units.length, real=True),
            Symbol("z", units.length, real=True),
        )

    def cartesian_transform(self, base_scalars: Optional[Sequence[Expr]] = None) -> Sequence[Expr]:
        return base_scalars or self.base_scalars

    def inverse_cartesian_transform(self, cartesian_scalars: Sequence[Expr]) -> Sequence[Expr]:
        return cartesian_scalars

    def generate_lame_coefficients(self) -> tuple[Expr, Expr, Expr]:
        return S.One, S.One, S.One

    def cartesian_derivative_matrix(
        self,
        _base_scalars: Optional[Sequence[Expr]] = None,
    ) -> Matrix:
        return Matrix.eye(3)

    def base_vector_matrix(
        self,
        _base_scalars: Optional[Sequence[Expr]] = None,
    ) -> Matrix:
        return Matrix.eye(3)

    def diff_base_vector_matrix(
        self,
        _wrt: SymSymbol,
        _base_scalars: Optional[Sequence[Expr]] = None,
    ) -> Matrix:
        return Matrix.zeros(3)


__all__ = ["CartesianCoordinateSystem"]

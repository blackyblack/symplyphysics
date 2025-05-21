from typing import Any, Optional

from sympy import Expr, Basic, Symbol as SymSymbol, S
from sympy.matrices.dense import DenseMatrix

from ..miscellaneous import sympify_expr
from ..vectors import is_vector_expr

from .base_system import BaseCoordinateSystem
from .point import AppliedPoint, check_point_with_system


class CoordinateScalar(Expr):
    _scalar: Expr
    _system: BaseCoordinateSystem
    _point: AppliedPoint | SymSymbol

    @property
    def scalar(self) -> Expr:
        return self._scalar

    @property
    def system(self) -> BaseCoordinateSystem:
        return self._system

    @property
    def point(self) -> AppliedPoint | SymSymbol:
        return self._point

    @property
    def args(self) -> tuple[Basic, ...]:
        return self.scalar, self.system, self.point

    def _hashable_content(self) -> tuple[Any, ...]:
        return self.args

    def __new__(
        cls,
        scalar: Any,
        system: BaseCoordinateSystem,
        point: Optional[AppliedPoint | SymSymbol] = None,
    ) -> Expr:
        if scalar == 0:
            return S.Zero

        scalar = sympify_expr(scalar)

        if isinstance(scalar, DenseMatrix):
            raise ValueError(f"Expected scalar, got matrix: {scalar}")

        if is_vector_expr(scalar):
            raise ValueError(f"Expected scalar, got vector: {scalar}")

        point = check_point_with_system(system, point)

        obj = super().__new__(cls)  # pylint: disable=no-value-for-parameter

        obj._scalar = scalar
        obj._system = system
        obj._point = point

        return obj


__all__ = ["CoordinateScalar"]

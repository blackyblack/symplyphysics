from typing import Any, Optional

from sympy import Expr, Basic, S
from sympy.matrices.dense import DenseMatrix
from sympy.physics.units import Dimension

from symplyphysics.core.dimensions.collect_quantity import collect_quantity_factor_and_dimension
from symplyphysics.core.symbols.symbols import BasicSymbol

from ..miscellaneous import sympify_expr
from ..vectors import is_vector_expr

from .base_system import BaseCoordinateSystem
from .point import AppliedPoint, check_point_with_system


class CoordinateScalar(Expr):
    """
    A `CoordinateScalar` represents the value of a scalar field at a particular *point* in space
    in a certain *coordinate system*. Common examples of scalar fields include the gravitational
    potential, the electric potential, temperature, pressure, etc.

    Note that unlike in `CoordinateVector`, the kind of coordinate system the scalar is defined in
    does not play a role in comparing them, i.e. both the system and the point of application must
    match.
    """

    _scalar: Expr
    _system: BaseCoordinateSystem
    _point: AppliedPoint | BasicSymbol

    @property
    def scalar(self) -> Expr:
        return self._scalar

    @property
    def system(self) -> BaseCoordinateSystem:
        return self._system

    @property
    def point(self) -> AppliedPoint | BasicSymbol:
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
        point: Optional[AppliedPoint | BasicSymbol] = None,
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

    def collect_quantity_factor_and_dimension(self) -> tuple[Expr, Dimension]:
        return collect_quantity_factor_and_dimension(self.scalar)


__all__ = ["CoordinateScalar"]

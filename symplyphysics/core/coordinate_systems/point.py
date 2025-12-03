from __future__ import annotations

from typing import Optional, Iterable, Any, Sized

from sympy import Basic, Expr, Dict, Tuple
from sympy.printing.printer import Printer

from symplyphysics.core.symbols.symbols import Symbol, BasicSymbol

from ..miscellaneous import sympify_expr

from .base_system import BaseCoordinateSystem
from .cartesian_system import CartesianCoordinateSystem


class AppliedPoint(Basic):
    """
    An `AppliedPoint` corresponds to a point in (3D) space whose coordinates are defined within a
    certain coordinate system.
    """

    _coordinates: Dict
    _system: BaseCoordinateSystem

    @property
    def coordinates(self) -> Dict:  # Dict[Symbol, Expr]
        return Dict(dict(zip(self.system.base_scalars, self.args[0])))

    @property
    def system(self) -> BaseCoordinateSystem:
        return self.args[1]

    def __getitem__(self, base_scalar: Symbol) -> Expr:
        return self.coordinates[base_scalar]

    def __new__(
        cls,
        coordinates: Iterable[Any],
        system: BaseCoordinateSystem,
    ) -> AppliedPoint:
        return super().__new__(cls)  # pylint: disable=no-value-for-parameter

    def __init__(
        self,
        coordinates: Iterable[Any],
        system: BaseCoordinateSystem,
    ):
        super().__init__()

        if isinstance(coordinates, Sized):
            n = len(coordinates)  # can't extract out of if-block, mypy complains otherwise
        else:
            coordinates = tuple(coordinates)
            n = len(coordinates)

        if n != 3:
            raise ValueError(f"The point must have all 3 coordinates defined, got {n}.")

        coordinates = Tuple(*(sympify_expr(coordinate) for coordinate in coordinates))

        self._args = coordinates, system

    def _sympystr(self, p: Printer) -> str:
        system_name = type(self.system).__name__.removesuffix("CoordinateSystem")
        point_name = f"{system_name}Point"

        coordinates = ", ".join(
            f"{p.doprint(s)} = {p.doprint(c)}" for s, c in self.coordinates.items())

        return f"{point_name}({coordinates})"

    def _hashable_content(self) -> tuple[Any, ...]:
        return self.args


# Used as a common point for Cartesian vectors
GLOBAL_POINT = BasicSymbol("P")


def check_point_with_system(
    system: BaseCoordinateSystem,
    point: Optional[AppliedPoint | BasicSymbol],
) -> AppliedPoint | BasicSymbol:
    if isinstance(system, CartesianCoordinateSystem):
        point = point or GLOBAL_POINT
    elif point is None:
        message = "The point of application must be defined for vectors in non-Cartesian systems"
        raise ValueError(message)

    if isinstance(point, AppliedPoint) and system != point.system:
        message = "The system of the vector and the point of its application must coincide"
        raise ValueError(message)

    return point


__all__ = [
    "AppliedPoint",
    "check_point_with_system",
    "GLOBAL_POINT",
]

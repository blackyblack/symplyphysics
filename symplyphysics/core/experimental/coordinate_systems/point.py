from __future__ import annotations

from typing import Mapping, Iterable, Any, Sized, SupportsFloat, Optional

from sympy import Basic, Expr, Symbol as SymSymbol
from sympy.printing.printer import Printer

from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.symbols.symbols import next_name

from ..miscellaneous import sympify_expr

from .base_system import BaseCoordinateSystem
from .cartesian_system import CartesianCoordinateSystem


def _prepare(
    coordinates: Iterable[SupportsFloat],
    system: BaseCoordinateSystem,
) -> dict[SymSymbol, Expr]:
    return {
        scalar: sympify_expr(coordinate)
        for scalar, coordinate in zip(system.base_scalars, coordinates)
    }


class AppliedPoint(Basic):
    """
    An `AppliedPoint` corresponds to a point in (3D) space whose coordinates are defined within a
    certain coordinate system.
    """

    _coordinates: Mapping[SymSymbol, Expr]
    _system: BaseCoordinateSystem

    # _iterable = False

    @property
    def coordinates(self) -> Mapping[SymSymbol, Expr]:
        return self._coordinates

    @property
    def system(self) -> BaseCoordinateSystem:
        return self._system

    def __getitem__(self, base_scalar: SymSymbol) -> Expr:
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

        self._coordinates = _prepare(coordinates, system)
        self._system = system

    def _sympystr(self, p: Printer) -> str:
        system_name = type(self.system).__name__.removesuffix("CoordinateSystem")
        point_name = f"{system_name}Point"

        coordinates = ", ".join(
            f"{p.doprint(s)} = {p.doprint(c)}" for s, c in self.coordinates.items())

        return f"{point_name}({coordinates})"

    def _hashable_content(self) -> tuple[Any, ...]:
        return (tuple(self.coordinates.items()),)

    def equals(self, other: AppliedPoint) -> bool:
        for base_scalar, coordinate in self.coordinates.items():
            if base_scalar not in other.coordinates:
                return False

            if not expr_equals(coordinate, other[base_scalar]):
                return False

        return True


# Used as a common point for Cartesian vectors
GLOBAL_POINT = SymSymbol(next_name("P"))


def check_point_with_system(
    system: BaseCoordinateSystem,
    point: Optional[AppliedPoint | SymSymbol],
) -> AppliedPoint | SymSymbol:
    if isinstance(system, CartesianCoordinateSystem):
        point = point or GLOBAL_POINT
    elif point is None:
        message = "The point of application must be defined for vectors in non-Cartesian systems"
        raise ValueError(message)

    if isinstance(point, AppliedPoint) and system != point.system:
        message = "The system of the vector and the point of its application must coincide"
        raise ValueError(message)

    return point


__all__ = ["AppliedPoint", "check_point_with_system"]

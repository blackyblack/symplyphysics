from __future__ import annotations

from collections.abc import Sized
from typing import Any, Iterable
from sympy import Expr, Symbol as SymSymbol, Atom
from sympy.printing.printer import Printer

from symplyphysics.core.expr_comparisons import expr_equals

from ..coordinate_systems.new_coordinate_systems import BaseCoordinateSystem
from ..miscellaneous import sympify_expr


def _prepare(coordinates: Iterable[Any], system: BaseCoordinateSystem) -> dict[SymSymbol, Expr]:
    return {
        scalar: sympify_expr(coordinate)
        for scalar, coordinate in zip(system.base_scalars, coordinates)
    }


class AppliedPoint(Atom):
    """
    An `AppliedPoint` corresponds to a point in (3D) space whose coordinates are defined within a
    certain coordinate system.
    """

    _coordinates: dict[SymSymbol, Expr]
    _system: BaseCoordinateSystem

    _iterable = False

    @property
    def coordinates(self) -> dict[SymSymbol, Expr]:
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


__all__ = [
    "AppliedPoint",
]

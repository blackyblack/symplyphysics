from __future__ import annotations

from typing import Optional, Any
from sympy import Basic, Expr, Atom, Symbol as SymSymbol, sympify
from sympy.printing.printer import Printer

from symplyphysics.core.symbols.id_generator import next_id

from ..coordinate_systems import BaseCoordinateSystem


class BasePoint(Basic):  # type: ignore[misc]
    """
    Base class for points. A point is an abstract geometrical entity that represents an exact
    position in physical (3D) space and that has no size.
    """


class PointSymbol(BasePoint, Atom):  # type: ignore[misc]
    """
    A `PointSymbol` corresponds to a certain point in (3D) space whose exact coordinates are not
    needed to be known.
    """

    _display_name: str
    _display_latex: str

    @property
    def display_name(self) -> str:
        return self._display_name

    @property
    def display_latex(self) -> str:
        return self._display_latex

    def __new__(
        cls,
        display_name: Optional[str] = None,
        *,
        display_latex: Optional[str] = None,
    ) -> PointSymbol:
        return Atom.__new__(cls)  # type: ignore[no-any-return]

    def __init__(
        self,
        display_name: Optional[str] = None,
        *,
        display_latex: Optional[str] = None,
    ) -> None:
        Atom.__init__(self)
        BasePoint.__init__(self)

        if display_name is None:
            id_ = next_id("PT")
            display_name = f"PT{id_}"
            if display_latex is None:
                display_latex = f"P_{{{id_}}}"
        elif display_latex is None:
            display_latex = display_name

        self._display_name = display_name
        self._display_latex = display_latex

    def _sympystr(self, _p: Printer) -> str:
        return self.display_name

    def _hashable_content(self) -> tuple[Any, ...]:
        return (id(self),)


def _prepare(coordinates: dict[SymSymbol, Any]) -> dict[SymSymbol, Expr]:
    return {scalar: sympify(coordinate, strict=True) for scalar, coordinate in coordinates.items()}


class AppliedPoint(BasePoint):
    """
    An `AppliedPoint` corresponds to a point in (3D) space whose coordinates are defined within a
    certain coordinate system.
    """

    _coordinates: dict[SymSymbol, Expr]
    _system: BaseCoordinateSystem

    @property
    def coordinates(self) -> dict[SymSymbol, Expr]:
        return self._coordinates

    @property
    def system(self) -> BaseCoordinateSystem:
        return self._system

    def __getitem__(self, scalar: SymSymbol) -> Expr:
        return self.coordinates[scalar]

    def __new__(
        cls,
        coordinates: dict[SymSymbol, Expr],
        system: BaseCoordinateSystem,
    ) -> AppliedPoint:
        return super().__new__(cls)  # type: ignore[no-any-return]

    def __init__(
        self,
        coordinates: dict[SymSymbol, Any],
        system: BaseCoordinateSystem,
    ):
        super().__init__()

        if set(coordinates.keys()) != set(system.base_scalars):
            raise ValueError("The point must have all coordinates defined.")

        self._coordinates = _prepare(coordinates)
        self._system = system

    def _sympystr(self, p: Printer) -> str:
        system_name = type(self.system).__name__.removesuffix("CoordinateSystem")
        point_name = f"{system_name}Point"

        coordinates = ", ".join(
            f"{p.doprint(s)} = {p.doprint(c)}" for s, c in self.coordinates.items())

        return f"{point_name}({coordinates})"

    def _hashable_content(self) -> tuple[Any, ...]:
        return (tuple(self.coordinates.items()),)


__all__ = [
    "BasePoint",
    "PointSymbol",
    "AppliedPoint",
]

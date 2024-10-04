"""
Points
======

A **point** in a 3D space is described by 3 coordinates which depend on the choice of the
coordinate system.
"""

from __future__ import annotations

from typing import Sequence, Iterator, Any

from sympy import Expr

from ...expr_comparisons import expr_equals
from ..coordinate_systems import BaseCoordinateSystem


class Point:
    """
    Class defining the `Point` functionality.
    """

    _coordinates: tuple[Expr, Expr, Expr]
    _system: BaseCoordinateSystem

    @property
    def coordinates(self) -> tuple[Expr, Expr, Expr]:
        """
        Coordinates of the point.
        """

        return self._coordinates

    @property
    def system(self) -> BaseCoordinateSystem:
        """
        Coordinate system the point is attached to.
        """

        return self._system

    def __init__(
        self,
        coordinates: Sequence[Expr],
        system: BaseCoordinateSystem,
    ) -> None:
        if len(coordinates) != 3:
            raise ValueError("Only 3-dimensional points are supported.")

        self._coordinates = coordinates[0], coordinates[1], coordinates[2]
        self._system = system

    def __getitem__(self, name: str) -> Expr:
        names = type(self.system).scalar_names()
        try:
            index = names.index(name)
        except ValueError as e:
            raise IndexError(
                f"Expected a scalar name from {names}, got {name!r}."
            ) from e

        return self.coordinates[index]

    def __iter__(self) -> Iterator[Expr]:
        return iter(self.coordinates)

    def __str__(self) -> str:
        sys_name = type(self.system).__name__
        index = sys_name.find("CoordinateSystem")
        name = f"{sys_name[:index]}Point"
        return f"{name}{self.coordinates}"

    def __eq__(self, other: object) -> bool:
        return (
            isinstance(other, Point)
            and self.system == other.system
            and all(expr_equals(x, y) for (x, y) in zip(self, other))
        )

    def subs(self, *args: Any) -> Point:
        """
        Substitute `old` for `new` in the point coordinates. See See `sympy.Basic.subs`
        for more info.
        """

        coordinates = [x.subs(*args) for x in self]
        return Point(coordinates, self.system)

    def simplify(self, **kwargs: Any) -> Point:
        """
        Simplify the point coordinates. See `sympy.simplify` for more info.
        """

        coordinates = [x.simplify(**kwargs) for x in self]
        return Point(coordinates, self.system)

    def evalf(self, **kwargs: Any) -> Point:
        """
        Evaluate the point coordinates to an accuracy of `kwargs["n"]` digits.
        """

        coordinates = [x.evalf(**kwargs) for x in self]
        return Point(coordinates, self.system)

    @staticmethod
    def from_system(system: BaseCoordinateSystem) -> Point:
        """
        Create a point whose coordinates are the given system's base scalars.
        """

        return Point(system.base_scalars, system)

    def apply_point(self, point: Point) -> Point:
        """
        Substitutes the base scalars of the given point's coordinate system with its
        coordinates within the current point's coordinates.
        """

        return self.subs(dict(zip(point.system.base_scalars, point)))

    @staticmethod
    def apply_to_expr(expr: Expr, point: Point) -> Expr:
        """
        Substitutes the base scalars of the given point's coordinate system with its
        coordinates within an expression.
        """

        return expr.subs(dict(zip(point.system.base_scalars, point)))

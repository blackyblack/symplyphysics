from typing import Iterable, TypeAlias
from sympy import Expr

Coordinate: TypeAlias = Expr | float


# This class not only represents point in space, but any trajectory in any-dimensional space.
# One should set it's coordinate to SymPy variable whenever it should be used as a trajectory, eg Point(C.x, C.y)
# represents a plane in 3D-space, where C is SympPy CoordSys3D.
class Point:
    # may contain not number but sympy expression, eg C.x
    _coordinates: list[Coordinate] = []

    def __init__(self, *coordinates: Coordinate) -> None:
        self._coordinates = list(coordinates)

    @property
    def coordinates(self) -> Iterable[Coordinate]:
        return iter(self._coordinates)

    def coordinate(self, index: int) -> Coordinate:
        if len(self._coordinates) <= index:
            return 0
        return self._coordinates[index]

    def set_coordinate(self, index: int, value: Coordinate) -> None:
        if value is None:
            value = 0
        if len(self._coordinates) <= index:
            self._coordinates.extend([0] * (index + 1 - len(self._coordinates)))
        self._coordinates[index] = value

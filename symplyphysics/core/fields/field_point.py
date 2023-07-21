from typing import Iterable, TypeAlias
from sympy import Expr

Coordinate: TypeAlias = Expr | float


# This class not only represents point in space, but any trajectory in any-dimensional space.
# One should set it's coordinate to SymPy variable whenever it should be used as a trajectory, eg FieldPoint(C.x, C.y)
# represents a plane in 3D-space, where C is SympPy CoordSys3D.
class FieldPoint:
    # may contain not number but sympy expression, eg C.x
    _coordinates: list[Coordinate] = []

    def __init__(self, *coordinates: Coordinate):
        self._coordinates = list(coordinates)

    @property
    def x(self) -> Coordinate:
        return self.coordinate(0)

    @x.setter
    def x(self, value_: Coordinate):
        self.set_coordinate(0, value_)

    @property
    def y(self) -> Coordinate:
        return self.coordinate(1)

    @y.setter
    def y(self, value_: Coordinate):
        self.set_coordinate(1, value_)

    @property
    def z(self) -> Coordinate:
        return self.coordinate(2)

    @z.setter
    def z(self, value_: Coordinate):
        self.set_coordinate(2, value_)

    @property
    def coordinates(self) -> Iterable[Coordinate]:
        return iter(self._coordinates)

    def coordinate(self, index: int) -> Coordinate:
        if len(self._coordinates) <= index:
            return 0
        return self._coordinates[index]

    def set_coordinate(self, index: int, value: Coordinate):
        if value is None:
            value = 0
        if len(self._coordinates) <= index:
            self._coordinates.extend([0] * (index + 1 - len(self._coordinates)))
        self._coordinates[index] = value

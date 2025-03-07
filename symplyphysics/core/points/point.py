from typing import Iterable, Any
from sympy import Expr, sympify, S


# This class not only represents point in space, but any trajectory in any-dimensional space.
# One should set it's coordinate to SymPy variable whenever it should be used as a trajectory, eg Point(C.x, C.y)
# represents a plane in 3D-space, where C is SympPy CoordSys3D.
class Point:
    # may contain not number but sympy expression, eg C.x
    _coordinates: list[Expr] = []

    def __init__(self, *coordinates: Any) -> None:
        self._coordinates = [sympify(c, strict=True) for c in coordinates]

    @property
    def coordinates(self) -> Iterable[Expr]:
        return iter(self._coordinates)

    def coordinate(self, index: int) -> Expr:
        if len(self._coordinates) <= index:
            return 0
        return self._coordinates[index]

    def set_coordinate(self, index: int, value: Any) -> None:
        if value is None:
            value = S.Zero
        if len(self._coordinates) <= index:
            self._coordinates.extend([S.Zero] * (index + 1 - len(self._coordinates)))
        self._coordinates[index] = sympify(value, strict=True)

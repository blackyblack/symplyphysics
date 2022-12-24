from typing import Any, List


# This class not only represents point in space, but any trajectory in any-dimensional space.
# One should set it's coordinate to SymPy variable whenever it should be used as a trajectory, eg FieldPoint(C.x, C.y)
# represents a plane in 3D-space, where C is SympPy CoordSys3D.
class FieldPoint:
    # may contain not number but sympy expression, eg C.x
    _coordinates: List[Any] = []

    def __init__(self, x_=0, y_=0, z_=0):
        self._coordinates = []
        if x_ != 0: self.set_coordinate(0, x_)
        if y_ != 0: self.set_coordinate(1, y_)
        if z_ != 0: self.set_coordinate(2, z_)

    @property
    def x(self):
        return self.coordinate(0)
    @property
    def y(self):
        return self.coordinate(1)
    @property
    def z(self):
        return self.coordinate(2)
    @x.setter
    def x(self, value_):
        self.set_coordinate(0, value_)
    @y.setter
    def y(self, value_):
        self.set_coordinate(1, value_)
    @z.setter
    def z(self, value_):
        self.set_coordinate(2, value_)

    @property
    def coordinates(self):
        return iter(self._coordinates)

    def coordinate(self, index: int):
        if len(self._coordinates) <= index: return 0
        return self._coordinates[index]

    def set_coordinate(self, index: int, value):
        if value is None: value = 0
        if len(self._coordinates) <= index:
            self._coordinates.extend([0] * (index + 1 - len(self._coordinates)))
        self._coordinates[index] = value

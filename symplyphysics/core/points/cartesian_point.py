from .point import Coordinate, Point


# This class represents point in 3d cartesian space (rectangular coordinates).
class CartesianPoint(Point):
    # Length of a rectangle
    @property
    def x(self) -> Coordinate:
        return self.coordinate(0)

    @x.setter
    def x(self, value_: Coordinate):
        self.set_coordinate(0, value_)

    # Width of a rectangle
    @property
    def y(self) -> Coordinate:
        return self.coordinate(1)

    @y.setter
    def y(self, value_: Coordinate):
        self.set_coordinate(1, value_)

    # Height of a rectangle
    @property
    def z(self) -> Coordinate:
        return self.coordinate(2)

    @z.setter
    def z(self, value_: Coordinate):
        self.set_coordinate(2, value_)

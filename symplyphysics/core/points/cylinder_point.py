from .point import Coordinate, Point


# This class represents point in 3d cylindrical space.
class CylinderPoint(Point):
    @property
    def radius(self) -> Coordinate:
        return self.coordinate(0)

    @radius.setter
    def radius(self, value_: Coordinate):
        self.set_coordinate(0, value_)

    @property
    def azimuthal_angle(self) -> Coordinate:
        return self.coordinate(1)

    @azimuthal_angle.setter
    def azimuthal_angle(self, value_: Coordinate):
        self.set_coordinate(1, value_)

    @property
    def height(self) -> Coordinate:
        return self.coordinate(2)

    @height.setter
    def height(self, value_: Coordinate):
        self.set_coordinate(2, value_)

    @property
    def r(self) -> Coordinate:
        return self.radius

    @r.setter
    def r(self, value_: Coordinate):
        self.radius = value_

    @property
    def theta(self) -> Coordinate:
        return self.azimuthal_angle

    @theta.setter
    def theta(self, value_: Coordinate):
        self.azimuthal_angle = value_

    @property
    def z(self) -> Coordinate:
        return self.height

    @z.setter
    def z(self, value_: Coordinate):
        self.height = value_

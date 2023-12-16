from .point import Coordinate, Point


# This class represents point in 3d spherical space.
class SpherePoint(Point):
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
    def polar_angle(self) -> Coordinate:
        return self.coordinate(2)

    @polar_angle.setter
    def polar_angle(self, value_: Coordinate):
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
    def phi(self) -> Coordinate:
        return self.polar_angle

    @phi.setter
    def phi(self, value_: Coordinate):
        self.polar_angle = value_

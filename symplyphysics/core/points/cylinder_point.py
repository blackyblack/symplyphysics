from typing import Any
from sympy import Expr
from .point import Point


# This class represents point in 3d cylindrical space.
class CylinderPoint(Point):

    @property
    def radius(self) -> Expr:
        return self.coordinate(0)

    @radius.setter
    def radius(self, value_: Any) -> None:
        self.set_coordinate(0, value_)

    @property
    def azimuthal_angle(self) -> Expr:
        return self.coordinate(1)

    @azimuthal_angle.setter
    def azimuthal_angle(self, value_: Any) -> None:
        self.set_coordinate(1, value_)

    @property
    def height(self) -> Expr:
        return self.coordinate(2)

    @height.setter
    def height(self, value_: Any) -> None:
        self.set_coordinate(2, value_)

    @property
    def r(self) -> Expr:
        return self.radius

    @r.setter
    def r(self, value_: Any) -> None:
        self.radius = value_

    @property
    def theta(self) -> Expr:
        return self.azimuthal_angle

    @theta.setter
    def theta(self, value_: Any) -> None:
        self.azimuthal_angle = value_

    @property
    def z(self) -> Expr:
        return self.height

    @z.setter
    def z(self, value_: Any) -> None:
        self.height = value_

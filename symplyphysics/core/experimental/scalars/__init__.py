"""
Scalars
=======

A **scalar** is a quantity that can be described by a single number or expression.
Scalars are unaffected by a coordinate system transformation and therefore do not
depend on them.

A **scalar field** is a mapping associating a single number to each point in space.
The ``Scalar`` class defined in this module provides a functionality to define both
a scalar and a scalar field.
"""

from __future__ import annotations

from typing import Optional, Any

from sympy import Expr

from ..points import Point


class Scalar:
    """
    Class defining the `Scalar` functionality.
    """

    _value: Expr
    _point: Optional[Point]

    def __init__(
        self,
        value: Expr,
        point: Optional[Point] = None,
    ) -> None:
        self._value = value
        self._point = point

    @property
    def value(self) -> Expr:
        """
        The value of the scalar.
        """

        return self._value

    @property
    def point(self) -> Optional[Point]:
        """
        The point in space the scalar describes. ``None`` if the scalar
        is not attached to any particular point.
        """

        return self._point

    def subs(self, *args: Any) -> Scalar:
        """
        Substitute `old` for `new` in the scalar value. See See `sympy.Basic.subs`
        for more info.
        """

        value = self.value.subs(*args)
        return Scalar(value, self.point)

    def simplify(self, **kwargs: Any) -> Scalar:
        """
        Simplify the scalar value. See `sympy.simplify` for more info.
        """

        value = self.value.simplify(**kwargs)
        return Scalar(value, self.point)

    def evalf(self, **kwargs: Any) -> Scalar:
        """
        Evaluate the scalar value to an accuracy of `kwargs["n"]` digits.
        """

        value = self.value.evalf(**kwargs)
        return Scalar(value, self.point)

    def apply_point(self, point: Point) -> Scalar:
        """
        Substitutes the base scalars of the given ``point``'s coordinate system with its
        coordinates within the scalar value..
        """

        return self.subs(dict(zip(point.system.base_scalars, point)))

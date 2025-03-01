"""
Vectors
=======

A **vector** is an object that has magnitude, or length, and direction. This module
introduces the functionality for constructing and operating with vectors.

A **vector field** is a mapping associating a vector to each point in a space, most
commonly the Euclidean space.
"""

from __future__ import annotations

from typing import Sequence, Optional, Any, Iterator

from sympy import Expr, S, sqrt, Matrix, zeros

from ..points import Point
from ..coordinate_systems import BaseCoordinateSystem, RectilinearCoordinateSystem


class Vector:
    """
    Class defining `Vector` functionality.
    """

    _components: Matrix
    _system: BaseCoordinateSystem
    _point: Optional[Point]

    def __init__(
        self,
        components: Sequence[Expr],
        system: BaseCoordinateSystem,
        point: Optional[Point] = None,
    ) -> None:
        """
        Create new vector.

        **Conditions:**

        - If the ``point`` is given, then it should be attached to the same coordinate
          system as the vector.
        """

        if point and point.system != system:
            raise ValueError(
                "The point of application must be in the "
                "same coordinate system as the vector."
            )

        if len(components) != 3:
            raise ValueError("Only 3-dimensional vectors are supported.")

        self._components = Matrix([components[0], components[1], components[2]])
        self._system = system
        self._point = point

    @property
    def components(self) -> Matrix:
        r"""
        Components of the vector with respect to the orthonormal basis:

        Formula:
            :code:`v = Sum(v[i] * b[i], i)`

        Latex:
            .. math::
                \mathbf{v} = \sum_i \hat{v}^i \mathbf{\hat{b}}_i

        **Notation:**

        - :math:`\mathbf{v}` (:code:`v`) is the vector in question.
        - :math:`\hat{v}^i` (:code:`v[i]`) is the i-th component of the vector.
        - :math:`\mathbf{\hat{b}}_i` (:code:`b[i]`) is the i-th normalized basis vector.
        """

        return self._components

    @property
    def system(self) -> BaseCoordinateSystem:
        """
        Coordinate system the vector is attached to.
        """

        return self._system

    @property
    def point(self) -> Optional[Point]:
        """
        The point the vector is applied to. ``None`` if there is no specific point of
        application.

        **Notes:**

        - The point is attached to the same coordinate system as the vector.
        """

        return self._point

    def __str__(self) -> str:
        sys_name = type(self.system).__name__
        index = sys_name.find("CoordinateSystem")
        name = f"{sys_name[:index]}Vector"
        v_1, v_2, v_3 = self
        return f"{name}{(v_1, v_2, v_3)}"

    def __getitem__(self, index: int) -> Expr:
        return self.components[index]

    def __iter__(self) -> Iterator[Expr]:
        return iter([self.components[0], self.components[1], self.components[2]])

    def __eq__(self, other: object) -> bool:
        return (
            isinstance(other, Vector)
            and self.system == other.system
            and self.point == other.point
            and (self - other).simplify().is_zero()
        )

    @staticmethod
    def zero(system: BaseCoordinateSystem, point: Optional[Point] = None) -> Vector:
        """
        Construct the zero vector. Note that the zero vector is independent of the coordinate system.
        """

        return Vector(zeros(3, 1), system, point)

    def is_zero(self) -> bool:
        """
        Check if the vector is zero.
        """

        return self.components.is_zero_matrix

    def assert_compatibility_with(self, other: Vector, *, strict: bool = True) -> None:
        """
        Raises ``ValueError`` iff ``self`` and ``other`` are either in different
        coordinate systems or applied at different points. If ``strict`` is set
        to ``False``, an exception is not raised if the vectors are applied at
        different points.
        """

        if self.system != other.system:
            raise ValueError("Vectors must be defined in the same coordinate system.")

        if self.is_zero() or other.is_zero():
            return

        if not strict and not isinstance(self.system, RectilinearCoordinateSystem):
            raise ValueError(
                "Vectors in curvilinear coordinate systems must be "
                "applied at the same point in order to be compatible."
            )

        if strict and self.point != other.point:
            raise ValueError("Vectors must be applied at the same point.")

    def add(self, other: Vector, *, strict: bool = True) -> Vector:
        """
        Adds two vectors. The vectors must be attached to the same coordinate system.
        They must also be applied at the same point unless ``strict`` is set to ``False``.
        """

        self.assert_compatibility_with(other, strict=strict)

        if self.is_zero():
            return other

        components = self.components + other.components
        return Vector(components, self.system, self.point)

    __add__ = add

    def __neg__(self) -> Vector:
        return S.NegativeOne * self

    def __pos__(self) -> Vector:
        return self

    def subtract(self, other: Vector, *, strict: bool = True) -> Vector:
        """
        Subtract two vectors. The vectors must be attached to the same coordinate system.
        They must also be applied at the same point unless ``strict`` is set to ``False``.
        """

        return self.add(-other, strict=strict)

    __sub__ = subtract

    def scale(self, scalar: Expr) -> Vector:
        """
        Scale the vector by a given expression.
        """

        components = scalar * self.components
        return Vector(components, self.system, self.point)

    __mul__ = scale
    __rmul__ = scale

    def __truediv__(self, scalar: Expr) -> Vector:
        return self * (S.One / scalar)

    def dot(self, other: Vector, *, strict: bool = True) -> Expr:
        """
        Perform the dot product between two vectors. The vectors must be attached to the same
        coordinate system. They must also be applied at the same point unless ``strict`` is
        set to ``False``.
        """

        self.assert_compatibility_with(other, strict=strict)

        return (self.components).dot(other.components)

    def cross(self, other: Vector, *, strict: bool = True) -> Vector:
        """
        Performs the cross product between two vectors. The vectors must be attached to the
        same coordinate system. They must also be applied at the same point unless ``strict``
        is set to ``False``.
        """

        self.assert_compatibility_with(other, strict=strict)

        return (self.components).cross(other.components)

    def magnitude(self) -> Expr:
        """
        Returns the magnitude of the vector.
        """

        return sqrt(self.dot(self))

    def __abs__(self) -> Vector:
        components = abs(self.components)
        return Vector(components, self.system, self.point)

    def normalize(self) -> Vector:
        """
        Performs vector normalization, i.e. divides the vector by its length, which yields
        a vector of unit length in the direction of the given vector.
        """

        return self / self.magnitude()

    def project_onto(self, target: Vector, *, strict: bool = True) -> Vector:
        """
        Projects ``self`` onto ``target``, i.e. extracts the component of ``self`` that is
        parallel to ``target``.
        """

        self_dot_target = self.dot(target, strict=strict)
        target_dot_target = target.dot(target, strict=strict)
        return (self_dot_target / target_dot_target) * target

    def reject_from(self, target: Vector, *, strict: bool = True) -> Vector:
        """
        Rejects ``self`` from ``target``, i.e. extracts the component of ``self`` that is
        orthogonal to ``target``.
        """

        return self - self.project_onto(target, strict=strict)

    def subs(self, *args: Any) -> Vector:
        r"""
        Substitute `old` for `new` in the vector components. See See `sympy.Basic.subs`
        for more info.
        """

        components = [a.subs(*args) for a in self]
        return Vector(components, self.system, self.point)

    def simplify(self, **kwargs: Any) -> Vector:
        r"""
        Simplify the vector components. See `sympy.simplify` for more info.
        """

        components = [a.simplify(**kwargs) for a in self]
        return Vector(components, self.system, self.point)

    def evalf(self, **kwargs: Any) -> Vector:
        """
        Evaluate the scalar value to an accuracy of `kwargs["n"]` digits.
        """

        components = [a.evalf(**kwargs) for a in self]
        return Vector(components, self.system, self.point)

    def apply_point(self, point: Point) -> Vector:
        """
        Substitutes the base scalars of the given ``point``'s coordinate system with its
        coordinates within the vector components.
        """

        return self.subs(dict(zip(point.system.base_scalars, point)))


__all__ = [
    "Vector",
]

from __future__ import annotations
from functools import partial
from typing import Callable, Sequence, TypeAlias, TypeVar
from sympy import Expr, sympify
from sympy.vector import express

from ..points.point import Point
from ..points.cartesian_point import CartesianPoint
from ..points.sphere_point import SpherePoint
from ..points.cylinder_point import CylinderPoint
from ..coordinate_systems.coordinate_systems import CoordinateSystem
from ...core.dimensions import ScalarValue

T = TypeVar("T", bound="Point")
FieldFunction: TypeAlias = Callable[[T], ScalarValue] | ScalarValue


def _subs_with_point(expr: ScalarValue, coordinate_system: CoordinateSystem,
    point_: Point) -> ScalarValue:
    base_scalars = coordinate_system.coord_system.base_scalars()
    # convert ScalarValue to Expr
    expression = sympify(expr)
    for i, scalar in enumerate(base_scalars):
        expression = expression.subs(scalar, point_.coordinate(i))
    return expression


# Contains mapping of point to a scalar value in _point_function, eg P(Point).
# Scalar field should not depend on the coordinate system, ie if applied to a point
# A1 in coordinate system C1 and having scalar value V as a result, it should
# have the same value V when applied to a point A2 in coordinate system C2, if
# A1 and A2 are the same points. Therefore scalar field is known to be invariant under
# Lorentz transformations.
class ScalarField:
    # Can contain lambda or some value. If value is stored, it will be returned when field is applied.
    _point_function: FieldFunction
    #NOTE: 4 and higher dimensional fields are not supported cause of using CoordSys3D
    #      to maintain ScalarField invariant.
    _coordinate_system: CoordinateSystem

    def __init__(self,
        point_function: FieldFunction = 0,
        coordinate_system: CoordinateSystem = CoordinateSystem(CoordinateSystem.System.CARTESIAN)):
        self._point_function = point_function
        self._coordinate_system = coordinate_system

    def __call__(self, point_: Point) -> ScalarValue:
        if not callable(self._point_function):
            return self._point_function
        # Point with general Point type is not checked against coordinate system.
        # It's up to user to make sure that field function works with general Point type.
        if isinstance(
                point_, CartesianPoint
        ) and self._coordinate_system.coord_system_type != CoordinateSystem.System.CARTESIAN:
            raise ValueError(
                f"Unsupported coordinate system for CartesianPoint: {self._coordinate_system}")
        if isinstance(
                point_, SpherePoint
        ) and self._coordinate_system.coord_system_type != CoordinateSystem.System.SPHERICAL:
            raise ValueError(
                f"Unsupported coordinate system for SpherePoint: {self._coordinate_system}")
        if isinstance(
                point_, CylinderPoint
        ) and self._coordinate_system.coord_system_type != CoordinateSystem.System.CYLINDRICAL:
            raise ValueError(
                f"Unsupported coordinate system for CylinderPoint: {self._coordinate_system}")
        return self._point_function(point_)

    @property
    def basis(self) -> Sequence[Expr]:
        return list(self.coordinate_system.coord_system.base_scalars())

    @property
    def coordinate_system(self) -> CoordinateSystem:
        return self._coordinate_system

    @property
    def field_function(self) -> FieldFunction:
        return self._point_function

    # Constructs new ScalarField from SymPy expression.
    # Can contain value instead of SymPy Vector, eg 0.5, but it should be sympified.
    @staticmethod
    def from_expression(
        expr: Expr | float,
        coordinate_system: CoordinateSystem = CoordinateSystem(CoordinateSystem.System.CARTESIAN)
    ) -> ScalarField:
        point_function = partial(_subs_with_point, expr, coordinate_system)
        return ScalarField(point_function, coordinate_system)

    # Applies field to a trajectory / surface / volume - calls field function with each element of the trajectory as parameter.
    # trajectory_ - list of expressions that correspond to a function in some space, eg [param, param] for a linear function y = x
    # return - value that depends on trajectory parameters.
    def apply(self, trajectory_: Sequence[Expr | float]) -> ScalarValue:
        trajectory_as_point = Point(*trajectory_)
        if self._coordinate_system.coord_system_type == CoordinateSystem.System.CARTESIAN:
            trajectory_as_point = CartesianPoint(*trajectory_)
        elif self._coordinate_system.coord_system_type == CoordinateSystem.System.SPHERICAL:
            trajectory_as_point = SpherePoint(*trajectory_)
        elif self._coordinate_system.coord_system_type == CoordinateSystem.System.CYLINDRICAL:
            trajectory_as_point = CylinderPoint(*trajectory_)
        return self(trajectory_as_point)

    # Convert coordinate system to space and apply field.
    # Applying field to entire space is necessary for SymPy field operators like Curl.
    # return - value that depends on basis parameters.
    def apply_to_basis(self) -> ScalarValue:
        return self.apply(self.basis)

    # Converts ScalarField to SymPy expression
    def to_expression(self) -> Expr:
        return sympify(self.apply_to_basis())

    # Convert field coordinate system to new basis and construct new field.
    # Scalar field invariant (coordinate system independence) should hold.
    def rebase(self, coordinate_system: CoordinateSystem) -> ScalarField:
        field_space_sympy = self.apply_to_basis()
        # Got a scalar value after applying to basis - use this value as field function
        if not isinstance(field_space_sympy, Expr):
            return ScalarField(field_space_sympy, coordinate_system)
        # Make linter happy
        field_space_expr: Expr = field_space_sympy
        if self.coordinate_system.coord_system_type != coordinate_system.coord_system_type:
            # This is a reverse transformation, if compared with Vector._extended_express()
            new_scalars = list(
                coordinate_system.transformation_to_system(
                self.coordinate_system.coord_system_type))
            for i, scalar in enumerate(self.coordinate_system.coord_system.base_scalars()):
                field_space_expr = field_space_expr.subs(scalar, new_scalars[i])
        # We do not want to maintain own field transformation functions, so
        # we convert our field to SymPy format, transform it and convert back to ScalarField.
        transformed_expr = express(field_space_expr,
            coordinate_system.coord_system,
            None,
            variables=True)
        return ScalarField.from_expression(transformed_expr, coordinate_system)

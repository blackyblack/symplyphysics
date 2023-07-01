from typing import Callable, Optional, TypeAlias
from sympy import Expr
from sympy.vector import express

from ..coordinate_systems.coordinate_systems import CoordinateSystem
from .field_point import FieldPoint

FieldFunction: TypeAlias = Callable[[FieldPoint], Expr] | Expr


# Converts SymPy expression to a lambda function to use in some field, eg
# converts expression C.x + C.y to lambda p: p.x + p.y
# If coordinate_system is not set, argument is returned as is, without conversion to lambda.
# Can contain value instead of SymPy Vector, eg 0.5.
def sympy_expression_to_field_function(
        expr: Expr, coordinate_system: Optional[CoordinateSystem] = None) -> FieldFunction:
    def _(point_: FieldPoint):
        # Duplicate check to make analyzer happy
        if coordinate_system is None:
            return expr
        base_scalars = coordinate_system.coord_system.base_scalars()
        # make a copy of expression
        expression = expr
        for i, scalar in enumerate(base_scalars):
            expression = expression.subs(scalar, point_.coordinate(i))
        return expression

    if coordinate_system is None:
        return expr
    return expr if not isinstance(expr, Expr) else _


# Contains mapping of point to a scalar value in _function, eg P(FieldPoint).
# Scalar field should not depend on the coordinate system, ie if applied to a point
# A1 in coordinate system C1 and having scalar value V as a result, it should
# have the same value V when applied to a point A2 in coordinate system C2, if
# A1 and A2 are the same points.
class ScalarField:
    # Can contain lambda or some value. If value is stored, it will be returned when field is applied.
    _point_function: FieldFunction
    #NOTE: 4 and higher dimensional fields are not supported cause of using CoordSys3D
    #      to maintain ScalarField invariant.
    _coordinate_system: Optional[CoordinateSystem] = None

    def __init__(self,
        point_function: FieldFunction = 0,
        coordinate_system: Optional[CoordinateSystem] = None):
        self._point_function = point_function
        self._coordinate_system = coordinate_system

    def __call__(self, point_: FieldPoint) -> Expr:
        return self._point_function(point_) if callable(
            self._point_function) else self._point_function

    @property
    def basis(self) -> list[Expr]:
        return list(self.coordinate_system.coord_system.base_scalars()
                   ) if self.coordinate_system is not None else []

    @property
    def coordinate_system(self) -> Optional[CoordinateSystem]:
        return self._coordinate_system

    @property
    def components(self) -> list[FieldFunction]:
        return [self._point_function]

    # Applies field to a trajectory / surface / volume - calls field function with each element of the trajectory as parameter.
    # trajectory_ - list of expressions that correspond to a function in some space, eg [param, param] for a linear function y = x
    # return - value that depends on trajectory parameters.
    def apply(self, trajectory_: list[Expr]) -> Expr:
        field_point = FieldPoint()
        for idx, element in enumerate(trajectory_):
            field_point.set_coordinate(idx, element)
        return self(field_point)

    # Convert coordinate system to space and apply field.
    # Applying field to entire space is necessary for SymPy field operators like Curl.
    # return - value that depends on basis parameters.
    def apply_to_basis(self) -> Expr:
        return self.apply(self.basis)


# Convert field coordinate system to new basis and construct new field.
# Scalar field invariant (coordinate system independence) should hold.
def field_rebase(field_: ScalarField,
    coordinate_system: Optional[CoordinateSystem] = None) -> ScalarField:
    # Simply set new coordinate system if field cannot be rebased
    if coordinate_system is None or field_.coordinate_system is None:
        field_function = 0 if len(field_.components) == 0 else field_.components[0]
        return ScalarField(field_function, coordinate_system)
    if coordinate_system.coord_system is None or field_.coordinate_system.coord_system is None:
        field_function = 0 if len(field_.components) == 0 else field_.components[0]
        return ScalarField(field_function, coordinate_system)
    return _extended_express(field_, coordinate_system)


def _extended_express(field_: ScalarField, system_to: CoordinateSystem) -> ScalarField:
    field_space_sympy = field_.apply_to_basis()
    if field_.coordinate_system.coord_system_type != system_to.coord_system_type:
        # This is a reverse transformation, if compared with Vector._extended_express()
        new_scalars = list(
            system_to.transformation_to_system(field_.coordinate_system.coord_system_type))
        for i, scalar in enumerate(field_.coordinate_system.coord_system.base_scalars()):
            field_space_sympy = field_space_sympy.subs(scalar, new_scalars[i])
    # We do not want to maintain own field transformation functions, so
    # we convert our field to SymPy format, transform it and convert back to ScalarField.
    transformed_vector_sympy = express(field_space_sympy,
        system_to.coord_system,
        None,
        variables=True)
    return field_from_sympy_vector(transformed_vector_sympy, system_to)


# Helpers to support SymPy field manipulations


# Constructs new ScalarField from SymPy expression using 'sympy_expression_to_field_function'.
def field_from_sympy_vector(sympy_vector_: Expr,
    coordinate_system: Optional[CoordinateSystem] = None) -> ScalarField:
    return ScalarField(sympy_expression_to_field_function(sympy_vector_, coordinate_system),
        coordinate_system)


# sympy_vector_from_field is identical to field.apply_to_basis()

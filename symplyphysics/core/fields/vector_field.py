from typing import Callable, Sequence, TypeAlias
from sympy import Expr, sympify
from sympy.vector import Vector as SymVector, express

from ..coordinate_systems.coordinate_systems import CoordinateSystem
from .field_point import FieldPoint
from ..vectors.vectors import Vector, vector_from_sympy_vector, sympy_vector_from_vector
from ...core.fields.scalar_field import ScalarValue

FieldFunction: TypeAlias = Callable[[FieldPoint], Sequence[ScalarValue]] | Sequence[ScalarValue]


# Converts SymPy expression to a lambda function to use in some field, eg
# converts expression C.x + C.y to lambda p: p.x + p.y
# If coordinate_system is not set, argument is returned as is, without conversion to lambda.
# Can contain value instead of SymPy Vector, eg 0.5.
def sympy_expression_to_field_function(expr: Sequence[ScalarValue],
    coordinate_system: CoordinateSystem) -> FieldFunction:

    def subs_with_point(point_: FieldPoint) -> Sequence[Expr]:
        # Duplicate check to make analyzer happy
        if coordinate_system is None:
            return expr
        base_scalars = coordinate_system.coord_system.base_scalars()
        result: list[Expr] = []
        for e in expr:
            expression = sympify(e)
            for i, scalar in enumerate(base_scalars):
                expression = expression.subs(scalar, point_.coordinate(i))
            result.append(expression)
        return result

    if coordinate_system is None:
        return expr
    return subs_with_point


# Contains mapping of point to vector in _point_function, eg P(FieldPoint).
# Vector field is coordinate system dependent, because generally vectors are not coordinate
# system invariant.
class VectorField:
    # Can contain lambda or some value. If value is stored, it will be returned when field is applied.
    _point_function: FieldFunction
    #NOTE: 4 and higher dimensional fields are not supported cause of using CoordSys3D
    #      that allows rebasing vector field to different coordinate systems.
    _coordinate_system: CoordinateSystem

    def __init__(self, coordinate_system: CoordinateSystem, point_function: FieldFunction):
        self._point_function = point_function
        self._coordinate_system = coordinate_system

    def __call__(self, point_: FieldPoint) -> Vector:
        result = self._point_function(point_) if callable(
            self._point_function) else self._point_function
        return Vector(self._coordinate_system, result)

    @property
    def basis(self) -> list[Expr]:
        return list(self._coordinate_system.coord_system.base_scalars())

    @property
    def coordinate_system(self) -> CoordinateSystem:
        return self._coordinate_system

    @property
    def field_function(self) -> FieldFunction:
        return self._point_function

    # Applies field to a trajectory / surface / volume - calls field functions with each element of the trajectory as parameter.
    # trajectory_ - list of expressions that correspond to a function in some space, eg [param, param] for a linear function y = x
    # return - vector parametrized by trajectory parameters.
    def apply(self, trajectory_: Sequence[Expr]) -> Vector:
        field_point = FieldPoint()
        for idx, element in enumerate(trajectory_):
            field_point.set_coordinate(idx, element)
        return self(field_point)

    # Convert coordinate system to space and apply field.
    # Applying field to entire space is necessary for SymPy field operators like Curl.
    # return - vector that depends on basis parameters.
    def apply_to_basis(self) -> Vector:
        return self.apply(self.basis)


# Convert field coordinate system to new basis and construct new field.
def field_rebase(field_: VectorField, coordinate_system: CoordinateSystem) -> VectorField:
    # Simply set new coordinate system if field cannot be rebased
    if (coordinate_system.coord_system is None or field_.coordinate_system.coord_system is None):
        return VectorField(coordinate_system, field_.field_function)

    field_space_sympy = sympy_vector_from_field(field_)
    if field_.coordinate_system.coord_system_type != coordinate_system.coord_system_type:
        # This is ScalarField._extended_express() but without transformation_to_system()
        new_scalars = list(coordinate_system.coord_system.base_scalars())
        for i, scalar in enumerate(field_.coordinate_system.coord_system.base_scalars()):
            field_space_sympy = field_space_sympy.subs(scalar, new_scalars[i])
    # We do not want to maintain own field transformation functions, so
    # we convert our field to SymPy format, transform it and convert back to VectorField.
    transformed_vector_sympy = express(field_space_sympy,
        coordinate_system.coord_system,
        None,
        variables=True)
    return field_from_sympy_vector(transformed_vector_sympy, coordinate_system)


# Helpers to support SymPy field manipulations, eg Curl


# Constructs new VectorField from SymPy expression using 'vector_from_sympy_vector'.
def field_from_sympy_vector(sympy_vector_: SymVector,
    coordinate_system: CoordinateSystem) -> VectorField:
    field_vector = vector_from_sympy_vector(sympy_vector_, coordinate_system)
    field_function = sympy_expression_to_field_function(field_vector.components, coordinate_system)
    return VectorField(coordinate_system, field_function)


# Apply field to entire coordinate system and convert to SymPy vector
def sympy_vector_from_field(field_: VectorField) -> SymVector:
    field_space = field_.apply_to_basis()
    return sympy_vector_from_vector(field_space)

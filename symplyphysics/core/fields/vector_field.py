from typing import Sequence
from sympy import Expr
from sympy.vector import Vector as SymVector, express

from ..coordinate_systems.coordinate_systems import CoordinateSystem
from .field_point import FieldPoint
from .scalar_field import FieldFunction, sympy_expression_to_field_function
from ..vectors.vectors import Vector, VectorComponent, vector_from_sympy_vector, sympy_vector_from_vector


# Contains mappings of point to vectors in components[], eg [P(FieldPoint), Q(FieldPoint)].
# Vector field can also be represented as a mapping of point to vector of values, eg
# lambda point: [point.x, point.y, point.z]. This way it is very similar to scalar field,
# where point is mapped to a value, eg lambda point: point.x
# Vector field should not depend on the coordinate system, ie if applied to a point
# A1 in coordinate system C1 and having vector V as a result, it should
# represent the same vector V when applied to a point A2 in coordinate system C2, if
# A1 and A2 are the same points.
class VectorField:
    #NOTE: 4 and higher dimensional fields are not supported cause of using CoordSys3D
    #      that allows rebasing vector field to different coordinate systems.
    _coordinate_system: CoordinateSystem
    # Can contain lambda or some value. If value is stored, it will be returned when field is applied.
    _components: list[FieldFunction] = []

    def __init__(self,
        coordinate_system: CoordinateSystem,
        vector_function_x_: FieldFunction = 0,
        vector_function_y_: FieldFunction = 0,
        vector_function_z_: FieldFunction = 0):
        self._coordinate_system = coordinate_system
        self._components = []
        if vector_function_x_ != 0:
            self.set_component(0, vector_function_x_)
        if vector_function_y_ != 0:
            self.set_component(1, vector_function_y_)
        if vector_function_z_ != 0:
            self.set_component(2, vector_function_z_)

    def __call__(self, point_: FieldPoint) -> Vector:
        vector_components: list[VectorComponent] = []
        for vector_function in self.components:
            vector_components.append(
                vector_function(point_) if callable(vector_function) else vector_function)
        return Vector(self._coordinate_system, vector_components)

    @property
    def basis(self) -> list[Expr]:
        return list(self._coordinate_system.coord_system.base_scalars()
                   ) if self._coordinate_system is not None else []

    @property
    def coordinate_system(self) -> CoordinateSystem:
        return self._coordinate_system

    @property
    def components(self) -> list[FieldFunction]:
        return self._components

    def component(self, index: int) -> FieldFunction:
        if len(self._components) <= index:
            return 0
        return self._components[index]

    # Useful for adding higher than 3 dimension field components
    def set_component(self, index: int, value: FieldFunction):
        if len(self._components) <= index:
            self._components.extend([0] * (index + 1 - len(self._components)))
        self._components[index] = value

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


# Constructs new VectorField from list and coordinate system
def _field_from_list(field_components_: Sequence[VectorComponent],
    coordinate_system: CoordinateSystem) -> VectorField:
    vector_components = [
        sympy_expression_to_field_function(c, coordinate_system) for c in field_components_
    ]
    if len(vector_components) < 3:
        vector_components.extend([0] * (3 - len(vector_components)))
    return VectorField(coordinate_system, vector_components[0], vector_components[1],
        vector_components[2])


# Convert field coordinate system to new basis and construct new field.
def field_rebase(field_: VectorField, coordinate_system: CoordinateSystem) -> VectorField:
    # Simply set new coordinate system if field cannot be rebased
    if (coordinate_system.coord_system is None or field_.coordinate_system.coord_system is None):
        vector_components = field_.components
        if len(vector_components) < 3:
            vector_components.extend([0] * (3 - len(vector_components)))
        return VectorField(coordinate_system, vector_components[0], vector_components[1],
            vector_components[2])

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
    return _field_from_list(field_vector.components, coordinate_system)


# Apply field to entire coordinate system and convert to SymPy vector
def sympy_vector_from_field(field_: VectorField) -> SymVector:
    field_space = field_.apply_to_basis()
    return sympy_vector_from_vector(field_space)

from typing import Any, List
from sympy.vector import Vector, CoordSys3D, express
from .field_point import FieldPoint
from .scalar_field import sympy_expression_to_field_function, extract_coord_system_from_sympy_vector
from ..vectors import array_to_sympy_vector, sympy_vector_to_array


# Contains mappings of point to vectors in components[], eg [P(FieldPoint), Q(FieldPoint)].
# Vector field can also be represented as a mapping of point to vector of values, eg
# lambda point: [point.x, point.y, point.z]. This way it is very similar to scalar field,
# where point is mapped to a value, eg lambda point: point.x
# Vector field should not depend on the coordinate system, ie if applied to a point
# A1 in coordinate system C1 and having vector V as a result, it should
# represent the same vector V when applied to a point A2 in coordinate system C2, if
# A1 and A2 are the same points.
class VectorField:
    # Can contain lambda or some value. If value is stored, it will be returned when field is applied.
    _components: List[Any] = []
    #NOTE: 4 and higher dimensional fields are not supported cause of using CoordSys3D
    #      that allows rebasing vector field to different coordinate systems.
    _coord_system: CoordSys3D = None

    def __init__(self, vector_function_x_=0, vector_function_y_=0, vector_function_z_=0, coord_system: CoordSys3D=None):
        self._components = []
        self._coord_system = coord_system
        if vector_function_x_ != 0: self.set_component(0, vector_function_x_)
        if vector_function_y_ != 0: self.set_component(1, vector_function_y_)
        if vector_function_z_ != 0: self.set_component(2, vector_function_z_)

    def __call__(self, point_: FieldPoint):
        result_vector = []
        for vector_function in self.components:
            result_vector.append(vector_function(point_) if callable(vector_function) else vector_function)
        return result_vector

    @property
    def basis(self) -> List[Any]:
        return list(self._coord_system.base_scalars()) if self._coord_system is not None else []

    @property
    def coord_system(self) -> CoordSys3D:
        return self._coord_system

    @property
    def components(self):
        return iter(self._components)

    def component(self, index: int) -> Any:
        if len(self._components) <= index: return 0
        return self._components[index]

    def set_component(self, index: int, value: Any):
        if len(self._components) <= index:
            self._components.extend([0] * (index + 1 - len(self._components)))
        self._components[index] = value

    # Applies field to a trajectory / surface / volume - calls field functions with each element of the trajectory as parameter.
    # trajectory_ - list of expressions that correspond to a function in some space, eg [param, param] for a linear function y = x
    # return - list that represents vector parametrized by trajectory parameters.
    def apply(self, trajectory_: List) -> List:
        field_point = FieldPoint()
        for idx, element in enumerate(trajectory_):
            field_point.set_coordinate(idx, element)
        return self(field_point)

    # Convert coordinate system to space and apply field.
    # Applying field to entire space is necessary for SymPy field operators like Curl.
    # return - vector that depends on basis parameters.
    def apply_to_basis(self) -> List[Any]:
        return self.apply(self.basis)

# Constructs new VectorField from SymPy expression using 'sympy_expression_to_field_function'.
def field_from_sympy_vector(sympy_vector_):
    field = VectorField()
    field_array = sympy_vector_to_array(sympy_vector_)
    field._coord_system = extract_coord_system_from_sympy_vector(sympy_vector_)
    for i in range(len(field_array)):
        field.set_component(i, sympy_expression_to_field_function(field_array[i]))
    return field

# Apply field to entire coordinate system and convert to SymPy vector
def field_to_sympy_vector(field_: VectorField) -> Vector:
    field_space = field_.apply_to_basis()
    return array_to_sympy_vector(field_.coord_system, field_space)

# Convert field coordinate system to new basis and construct new field.
def field_rebase(field_: VectorField, coord_system: CoordSys3D=None) -> VectorField:
    # Simply set new coordinate system if field cannot be rebased
    if coord_system is None or field_.coord_system is None:
        field = VectorField(0, 0, 0, coord_system)
        for idx, c in enumerate(field_.components):
            field.set_component(idx, c)
        return field
    # We do not want to maintain own vector transformation functions, so
    # we convert our field to SymPy format, transform it and convert back to VectorField.
    field_space_sympy = field_to_sympy_vector(field_)
    transformed_field_space = express(field_space_sympy, coord_system, variables=True)
    return field_from_sympy_vector(transformed_field_space)

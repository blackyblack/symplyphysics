from types import FunctionType
from typing import List
from sympy.vector import CoordSys3D, Vector
from .field_point import FieldPoint
from .scalar_field import sympy_expression_to_field_function, extract_coord_system_from_sympy_vector
from ..vectors import sympy_vector_to_array


# Contains mappings of point to vectors in components[], eg [P(FieldPoint), Q(FieldPoint)].
# Therefore vector field can be represented as vector of scalar fields.
# Vector field can also be represented as a mapping of point to vector of values, eg
# lambda point: [point.x, point.y, point.z]. This way it is very similar to scalar field,
# where point is mapped to a value, eg lambda point: point.x
class VectorField:
    _components: List[FunctionType] = []

    def __init__(self, vector_function_x_=0, vector_function_y_=0, vector_function_z_=0):
        self._components = []
        if vector_function_x_ != 0: self.set_component(0, vector_function_x_)
        if vector_function_y_ != 0: self.set_component(1, vector_function_y_)
        if vector_function_z_ != 0: self.set_component(2, vector_function_z_)

    @property
    def components(self):
        return iter(self._components)

    def component(self, index: int):
        if len(self._components) <= index: return 0
        return self._components[index]

    def set_component(self, index: int, value):
        if len(self._components) <= index:
            self._components.extend([0] * (index + 1 - len(self._components)))
        self._components[index] = value


# Converts SymPy vector, eg C.x * C.i + C.y * C.j, where C is CoordSys3D, to a field
def sympy_vector_to_field(sympy_vector_: Vector) -> VectorField:
    field = VectorField()
    field_array = sympy_vector_to_array(sympy_vector_)
    coord_system = extract_coord_system_from_sympy_vector(sympy_vector_)
    for i in range(len(field_array)):
        field.set_component(i, sympy_expression_to_field_function(coord_system, field_array[i]))
    return field

# Applies field to a trajectory / surface / volume - calls field functions with each element of the trajectory as parameter.
# field_ - VectorField with functions that map base coordinates to a vector. Can contain 0 instead of empty function.
# trajectory_ - list of expressions that correspond to a function in some space, eg [param, param] for a linear function y = x
# return - list that represents vector parametrized by trajectory parameters.
def apply_field(field_: VectorField, trajectory_: List) -> List:
    field_point = FieldPoint()
    for idx, element in enumerate(trajectory_):
        field_point.set_coordinate(idx, element)
    result_vector = []
    for vector_function in field_.components:
        result_vector.append(vector_function(field_point) if callable(vector_function) else vector_function)
    return result_vector

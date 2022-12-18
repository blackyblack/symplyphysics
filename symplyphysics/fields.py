from types import FunctionType
from typing import Any, List
from sympy import Expr
from sympy.vector import CoordSys3D, Vector
from sympy.vector.operators import _get_coord_systems

from symplyphysics.vectors import sympy_vector_to_array


class FieldPoint:
    # may contain not number but sympy expression, eg C.x
    _coordinates: List[Any] = []

    def __init__(self, x_=0, y_=0, z_=0):
        self._coordinates = []
        if x_ != 0: self.set_coordinate(0, x_)
        if y_ != 0: self.set_coordinate(1, y_)
        if z_ != 0: self.set_coordinate(2, z_)

    @property
    def x(self):
        return self.coordinate(0)
    @property
    def y(self):
        return self.coordinate(1)
    @property
    def z(self):
        return self.coordinate(2)
    @x.setter
    def x(self, value_):
        self.set_coordinate(0, value_)
    @y.setter
    def y(self, value_):
        self.set_coordinate(1, value_)
    @z.setter
    def z(self, value_):
        self.set_coordinate(2, value_)

    def coordinate(self, index: int):
        if len(self._coordinates) <= index: return 0
        return self._coordinates[index]

    def set_coordinate(self, index: int, value):
        if len(self._coordinates) <= index:
            self._coordinates.extend([0] * (index + 1 - len(self._coordinates)))
        self._coordinates[index] = value


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
    # Better replacement for lambda
    def subs_with_point(coord_system_: CoordSys3D, expr_):
        def _(point_: FieldPoint):
            base_scalars = coord_system_.base_scalars()
            # make a copy of expression
            expr = expr_
            for i in range(len(base_scalars)):
                expr = expr.subs(base_scalars[i], point_.coordinate(i))
            return expr

        if not isinstance(expr_, Expr):
            return expr_
        return expr_ if coord_system_ is None else _

    field = VectorField()
    field_array = sympy_vector_to_array(sympy_vector_)
    coord_system_set = _get_coord_systems(sympy_vector_)
    coord_system = None if len(coord_system_set) == 0 else next(iter(coord_system_set))
    for i in range(len(field_array)):
        field.set_component(i, subs_with_point(coord_system, field_array[i]))
    return field

# Convert coordinate system to space. Field can be applied to this space, eg apply_field(field, space).
# Applying field to entire space is necessary for SymPy field operators like Curl.
def coord_system_to_space(coord_system_: CoordSys3D) -> List:
    return list(coord_system_.base_scalars())

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

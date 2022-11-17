from types import FunctionType
from typing import Any, List
from sympy.vector import CoordSys3D


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


# Contains mappings of point to vectors in components[], eg [P(FieldPoint), Q(FieldPoint)]
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


    # Applies field to a trajectory (surface) - calls field functions with each element of the trajectory as parameter.
    # coord_system_ - CoordSys3D coordinate system to work with. Supports Cartesian, Spherical, Cylindrical coordinates.
    # field_ - VectorField with functions that map base coordinates to a vector. Each function should have the same number of 
    #          parameters as base coordinates in coord_system_. Can contain 0 instead of empty function.
    # trajectory_ - list of expressions that correspond to a function in some space. Each element of list
    #               corresponds to base coordinate in coord_system_.
    # return - list that represents vector parametrized by trajectory parameters.
    def apply(self, coord_system_: CoordSys3D, trajectory_: List) -> List:
        base_scalars = coord_system_.base_scalars()
        dimensions = len(base_scalars)
        field_point = FieldPoint()
        for i in range(min(dimensions, len(trajectory_))):
            field_point.set_coordinate(i, trajectory_[i])
        result_vector = []
        for vector_function in self._components:
            result_vector.append(vector_function(field_point) if callable(vector_function) else vector_function)
        return result_vector


# Converts vector with unit scalars (C.x, C.y) as parameters to field
def field_from_unit_vector(coord_system_: CoordSys3D, vector_: List) -> VectorField:
    # Better replacement for lambda
    def subs_with_point(coord_system_: CoordSys3D, expr_):
        def _(point_: FieldPoint):
            base_scalars = coord_system_.base_scalars()
            # make a copy of expression
            expr = expr_
            for i in range(len(base_scalars)):
                expr = expr.subs(base_scalars[i], point_.coordinate(i))
            return expr
        return _

    base_scalars = coord_system_.base_scalars()
    dimensions = len(base_scalars)
    field = VectorField()
    for i in range(min(dimensions, len(vector_))):
        field.set_component(i, subs_with_point(coord_system_, vector_[i]))
    return field

from types import FunctionType
from typing import Any, List
from sympy import Expr
from sympy.vector import CoordSys3D, Vector
from sympy.vector.operators import _get_coord_systems
from .field_point import FieldPoint


# Contains mapping of point to a scalar value in _function, eg P(FieldPoint).
class ScalarField:
    _function: FunctionType = None

    def __init__(self, scalar_function_=0):
        self._function = scalar_function_

    def __call__(self, point_: FieldPoint):
        return self._function(point_) if callable(self._function) else self._function


# Converts SymPy Vector, eg C.x + C.y to lambda with FieldPoint as parameter, eg
# lambda point: point.x + point.y
# Can contain value instead of SymPy Vector, eg 0.5.
def sympy_expression_to_field_function(coord_system_: CoordSys3D, expression_):
    def _(point_: FieldPoint):
        base_scalars = coord_system_.base_scalars()
        # make a copy of expression
        expression = expression_
        for i in range(len(base_scalars)):
            expression = expression.subs(base_scalars[i], point_.coordinate(i))
        return expression

    if not isinstance(expression_, Expr):
        return expression_
    return expression_ if coord_system_ is None else _

# Detects coordinate system being used in SymPy Vector.
# Throws an exception if multiple coordinate systems are detected.
def extract_coord_system_from_sympy_vector(sympy_vector_: Vector) -> CoordSys3D:
    coord_system_set = _get_coord_systems(sympy_vector_)
    if len(coord_system_set) > 1:
        coord_sys_names = [str(c) for c in coord_system_set]
        raise TypeError(f"Different coordinate systems in expression: {str(coord_sys_names)}")
    return None if len(coord_system_set) == 0 else next(iter(coord_system_set))

# Converts SymPy vector, eg C.x + C.y, where C is CoordSys3D, to a field
def sympy_vector_to_field(sympy_vector_: Vector) -> ScalarField:
    coord_system = extract_coord_system_from_sympy_vector(sympy_vector_)
    return ScalarField(sympy_expression_to_field_function(coord_system, sympy_vector_))

# Applies field to a trajectory / surface / volume - calls field function with each element of the trajectory as parameter.
# field_ - ScalarField with function that maps FieldPoint to a scalar value. Can contain 0 instead of empty function.
# trajectory_ - list of expressions that correspond to a function in some space, eg [param, param] for a linear function y = x
# return - value that depends on trajectory parameters.
def apply_field(field_: ScalarField, trajectory_: List) -> Any:
    field_point = FieldPoint()
    for idx, element in enumerate(trajectory_):
        field_point.set_coordinate(idx, element)
    return field_(field_point)

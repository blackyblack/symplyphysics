from functools import partial
from types import FunctionType
from typing import Any, List
from sympy import Expr
from sympy.vector import CoordSys3D, Vector, express
from sympy.vector.operators import _get_coord_systems
from .field_point import FieldPoint


# Detects coordinate system being used in SymPy Vector.
# Throws an exception if multiple coordinate systems are detected.
def extract_coord_system_from_sympy_vector(sympy_vector_: Vector) -> CoordSys3D:
    coord_system_set = _get_coord_systems(sympy_vector_)
    if len(coord_system_set) > 1:
        coord_sys_names = [str(c) for c in coord_system_set]
        raise TypeError(f"Different coordinate systems in expression: {str(coord_sys_names)}")
    return None if len(coord_system_set) == 0 else next(iter(coord_system_set))

def sympy_expression_to_field_function(sympy_vector_) -> Any:
    def _(coord_system_: CoordSys3D, point_: FieldPoint):
        base_scalars = coord_system_.base_scalars()
        # make a copy of expression
        expression = sympy_vector_
        for i in range(len(base_scalars)):
            expression = expression.subs(base_scalars[i], point_.coordinate(i))
        return expression

    coord_system = extract_coord_system_from_sympy_vector(sympy_vector_)
    if coord_system is None:
        return sympy_vector_
    return sympy_vector_ if not isinstance(sympy_vector_, Expr) else partial(_, coord_system)


# Contains mapping of point to a scalar value in _function, eg P(FieldPoint).
# Scalar field should not depend on the coordinate system, ie if applied to a point
# A1 in coordinate system C1 and having scalar value V as a result, it should
# have the same value V when applied to a point A2 in coordinate system C2, if
# A1 and A2 are the same points.
class ScalarField:
    _point_function: FunctionType = None
    #NOTE: 4 and higher dimensional fields are not supported cause of using CoordSys3D
    #      to maintain ScalarField invariant.
    _coord_system: CoordSys3D = None

    def __init__(self, point_function=0, coord_system: CoordSys3D=None):
        self._point_function = point_function
        self._coord_system = coord_system

    def __call__(self, point_: FieldPoint):
        return self._point_function(point_) if callable(self._point_function) else self._point_function

    @property
    def basis(self):
        return self._coord_system.base_scalars() if self._coord_system is not None else []

    # Converts SymPy Vector, eg C.x + C.y to a field with field function as a lambda with FieldPoint as parameter, eg
    # lambda point: point.x + point.y
    # Can contain value instead of SymPy Vector, eg 0.5.
    def from_sympy_vector(sympy_vector_):
        field = ScalarField()
        field._point_function = sympy_expression_to_field_function(sympy_vector_)
        field._coord_system = extract_coord_system_from_sympy_vector(sympy_vector_)
        return field

    # Applies field to a trajectory / surface / volume - calls field function with each element of the trajectory as parameter.
    # trajectory_ - list of expressions that correspond to a function in some space, eg [param, param] for a linear function y = x
    # return - value that depends on trajectory parameters.
    def apply(self, trajectory_: List) -> Any:
        field_point = FieldPoint()
        for idx, element in enumerate(trajectory_):
            field_point.set_coordinate(idx, element)
        return self._point_function(field_point)

    # Convert coordinate system to space and apply field.
    # Applying field to entire space is necessary for SymPy field operators like Curl.
    def apply_to_basis(self) -> List:
        return self.apply(list(self.basis))

    # Convert field coordinate system to new basis.
    # Scalar field invariant (coordinate system independence) should hold.
    def rebase_field(self, coord_system: CoordSys3D=None):
        if coord_system is None:
            raise TypeError(f"Rebasing to empty coordinate system is not supported")
        if self._coord_system is None:
            raise TypeError(f"Rebasing from empty coordinate system is not supported")
        field_space = self.apply_to_basis()
        transformed_field_space = express(field_space, coord_system, variables=True)
        transformed_field = ScalarField.from_sympy_vector(transformed_field_space)
        self.__init__(transformed_field._point_function, transformed_field._coord_system)

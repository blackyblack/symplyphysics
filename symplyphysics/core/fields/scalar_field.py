from typing import Any, List
from sympy import Expr
from sympy.vector import CoordSys3D, express

from symplyphysics.core.coordinate_systems.coordinate_systems import CoordinateSystem
from .field_point import FieldPoint


# Converts SymPy expression to a lambda function to use in some field, eg
# converts expression C.x + C.y to lambda p: p.x + p.y
# If coordinate_system is not set, argument is returned as is, without conversion to lambda.
# Can contain value instead of SymPy Vector, eg 0.5.
def sympy_expression_to_field_function(sympy_vector_: Any, coordinate_system: CoordinateSystem=None) -> Any:
    def _(point_: FieldPoint):
        base_scalars = coordinate_system.coord_system.base_scalars()
        # make a copy of expression
        expression = sympy_vector_
        for i in range(len(base_scalars)):
            expression = expression.subs(base_scalars[i], point_.coordinate(i))
        return expression

    if coordinate_system is None:
        return sympy_vector_
    return sympy_vector_ if not isinstance(sympy_vector_, Expr) else _


# Contains mapping of point to a scalar value in _function, eg P(FieldPoint).
# Scalar field should not depend on the coordinate system, ie if applied to a point
# A1 in coordinate system C1 and having scalar value V as a result, it should
# have the same value V when applied to a point A2 in coordinate system C2, if
# A1 and A2 are the same points.
class ScalarField:
    # Can contain lambda or some value. If value is stored, it will be returned when field is applied.
    _point_function: Any = None
    #NOTE: 4 and higher dimensional fields are not supported cause of using CoordSys3D
    #      to maintain ScalarField invariant.
    _coordinate_system: CoordinateSystem = None

    def __init__(self, point_function=0, coord_system: CoordSys3D=None):
        self._point_function = point_function
        self._coord_system = coord_system

    def __call__(self, point_: FieldPoint):
        return self._point_function(point_) if callable(self._point_function) else self._point_function

    @property
    def basis(self) -> List[Any]:
        return list(self.coordinate_system.coord_system.base_scalars()) if self.coordinate_system is not None else []

    @property
    def coordinate_system(self) -> CoordinateSystem:
        return self._coordinate_system

    @property
    def components(self):
        return [self._point_function]

    # Applies field to a trajectory / surface / volume - calls field function with each element of the trajectory as parameter.
    # trajectory_ - list of expressions that correspond to a function in some space, eg [param, param] for a linear function y = x
    # return - value that depends on trajectory parameters.
    def apply(self, trajectory_: List) -> Any:
        field_point = FieldPoint()
        for idx, element in enumerate(trajectory_):
            field_point.set_coordinate(idx, element)
        return self(field_point)

    # Convert coordinate system to space and apply field.
    # Applying field to entire space is necessary for SymPy field operators like Curl.
    # return - value that depends on basis parameters.
    def apply_to_basis(self) -> Any:
        return self.apply(self.basis)

# Constructs new ScalarField from SymPy expression using 'sympy_expression_to_field_function'.
def field_from_sympy_vector(sympy_vector_, coordinate_system: CoordinateSystem):
    return ScalarField(sympy_expression_to_field_function(sympy_vector_, coordinate_system), coordinate_system)

# Convert field coordinate system to new basis and construct new field.
# Scalar field invariant (coordinate system independence) should hold.
def field_rebase(field_: ScalarField, coordinate_system: CoordinateSystem=None) -> ScalarField:
    # Simply set new coordinate system if field cannot be rebased
    if coordinate_system is None or field_.coordinate_system is None:
        field_function = 0 if len(field_.components) == 0 else field_.components[0]
        return ScalarField(field_function, coordinate_system)

    # Similar to VectorField field_rebase but cannot implement generic rebase algorithm.
    field_space_sympy = field_.apply_to_basis()
    transformed_field_space = express(field_space_sympy, coordinate_system, variables=True)
    return field_from_sympy_vector(transformed_field_space)

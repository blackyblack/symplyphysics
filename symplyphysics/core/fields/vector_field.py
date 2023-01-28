from typing import Any, List
from sympy.vector import Vector as SympyVector, CoordSys3D
from .field_point import FieldPoint
from .scalar_field import sympy_expression_to_field_function
from ..vectors.vectors import Vector, extract_coord_system_from_sympy_vector, vector_from_sympy_vector, sympy_vector_from_vector, vector_rebase


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

    def __call__(self, point_: FieldPoint) -> Vector:
        vector_components = []
        for vector_function in self.components:
            vector_components.append(vector_function(point_) if callable(vector_function) else vector_function)
        return Vector(vector_components, self._coord_system)

    @property
    def basis(self) -> List[Any]:
        return list(self._coord_system.base_scalars()) if self._coord_system is not None else []

    @property
    def coord_system(self) -> CoordSys3D:
        return self._coord_system

    @property
    def components(self):
        return self._components

    def component(self, index: int) -> Any:
        if len(self._components) <= index: return 0
        return self._components[index]

    def set_component(self, index: int, value: Any):
        if len(self._components) <= index:
            self._components.extend([0] * (index + 1 - len(self._components)))
        self._components[index] = value

    # Applies field to a trajectory / surface / volume - calls field functions with each element of the trajectory as parameter.
    # trajectory_ - list of expressions that correspond to a function in some space, eg [param, param] for a linear function y = x
    # return - vector parametrized by trajectory parameters.
    def apply(self, trajectory_: List) -> Vector:
        field_point = FieldPoint()
        for idx, element in enumerate(trajectory_):
            if self._coord_system is not None:
                element_coord_system = extract_coord_system_from_sympy_vector(element)
                if element_coord_system is not None and element_coord_system != self._coord_system:
                    raise TypeError(f"Different coordinate systems in field and expression: {str(self._coord_system)} vs {str(element_coord_system)}")
            field_point.set_coordinate(idx, element)
        return self(field_point)

    # Convert coordinate system to space and apply field.
    # Applying field to entire space is necessary for SymPy field operators like Curl.
    # return - vector that depends on basis parameters.
    def apply_to_basis(self) -> Vector:
        return self.apply(self.basis)

# Constructs new VectorField from Vector
def field_from_vector(vector_: Vector) -> VectorField:
    vector_components = []
    for component in vector_.components:
        vector_components.append(sympy_expression_to_field_function(component, vector_.coord_system))
    if len(vector_components) < 3:
        vector_components.extend([0] * (3 - len(vector_components)))
    return VectorField(vector_components[0], vector_components[1], vector_components[2], vector_.coord_system)

# Convert field coordinate system to new basis and construct new field.
def field_rebase(field_: VectorField, coord_system: CoordSys3D=None) -> VectorField:
    # Simply set new coordinate system if field cannot be rebased
    if coord_system is None or field_.coord_system is None:
        return field_from_vector(Vector(field_.components, coord_system))
    # Transform field to vector and use vector_rebase to rebase field. Transform
    # back to field after rebasing.
    field_space = field_.apply_to_basis()
    transformed_field_space = vector_rebase(field_space, coord_system)
    return field_from_vector(transformed_field_space)

# Helpers to support SymPy field manipulations, eg Curl

# Constructs new VectorField from SymPy expression using 'vector_from_sympy_vector'.
def field_from_sympy_vector(sympy_vector_: Any) -> VectorField:
    field_vector = vector_from_sympy_vector(sympy_vector_)
    return field_from_vector(field_vector)

# Apply field to entire coordinate system and convert to SymPy vector
def sympy_vector_from_field(field_: VectorField) -> SympyVector:
    field_space = field_.apply_to_basis()
    return sympy_vector_from_vector(field_space)

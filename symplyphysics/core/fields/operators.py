from sympy import diff
from ..fields.vector_field import VectorField
from ..coordinate_systems.coordinate_systems import CoordinateSystem
from ..vectors.vectors import Vector


# Calculate Curl of the field, which is Cross(Nabla, Field)
def curl_operator(field_: VectorField) -> VectorField:
    if field_.coordinate_system.coord_system_type != CoordinateSystem.System.CARTESIAN:
        coord_name_from = CoordinateSystem.system_to_transformation_name(
            field_.coordinate_system.coord_system_type)
        raise ValueError(f"Curl is only supported for cartesian coordinates: got {coord_name_from}")
    field_space = field_.apply_to_basis()
    dimensions = 3
    if len(field_space.components) > dimensions:
        raise ValueError(
            f"Curl is only defined for {dimensions} dimensions. Got: {len(field_space.components)}")
    x = field_space.coordinate_system.coord_system.base_scalars()[0]
    y = field_space.coordinate_system.coord_system.base_scalars()[1]
    z = field_space.coordinate_system.coord_system.base_scalars()[2]
    field_x = field_space.components[0]
    field_y = field_space.components[1]
    field_z = field_space.components[2]
    field_rotor_vector = Vector([
        diff(field_z, y) - diff(field_y, z),
        diff(field_x, z) - diff(field_z, x),
        diff(field_y, x) - diff(field_x, y)
    ], field_space.coordinate_system)
    #TODO: avoid conversion to SymPy vector
    field_rotor_vector_sympy = field_rotor_vector.to_sympy_vector()
    return VectorField.from_sympy_vector(field_rotor_vector_sympy, field_space.coordinate_system)

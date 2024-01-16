from sympy import diff, sin, tan

from ..dimensions import ScalarValue
from ..fields.vector_field import VectorField
from ..fields.scalar_field import ScalarField
from ..coordinate_systems.coordinate_systems import CoordinateSystem
from ..vectors.vectors import Vector


def gradient_operator(field: ScalarField) -> Vector:
    field_space = field.apply_to_basis()
    if field.coordinate_system.coord_system_type == CoordinateSystem.System.CARTESIAN:
        x = field.coordinate_system.coord_system.base_scalars()[0]
        y = field.coordinate_system.coord_system.base_scalars()[1]
        z = field.coordinate_system.coord_system.base_scalars()[2]
        gradient = Vector([
            diff(field_space, x),
            diff(field_space, y),
            diff(field_space, z),
        ], field.coordinate_system)
        return gradient
    if field.coordinate_system.coord_system_type == CoordinateSystem.System.CYLINDRICAL:
        r = field.coordinate_system.coord_system.base_scalars()[0]
        theta = field.coordinate_system.coord_system.base_scalars()[1]
        z = field.coordinate_system.coord_system.base_scalars()[2]
        gradient = Vector([
            diff(field_space, r),
            diff(field_space, theta) / r,
            diff(field_space, z),
        ], field.coordinate_system)
        return gradient
    if field.coordinate_system.coord_system_type == CoordinateSystem.System.SPHERICAL:
        r = field.coordinate_system.coord_system.base_scalars()[0]
        phi = field.coordinate_system.coord_system.base_scalars()[1]
        theta = field.coordinate_system.coord_system.base_scalars()[2]
        gradient = Vector([
            diff(field_space, r),
            diff(field_space, phi) / (r * sin(theta)),
            diff(field_space, theta) / r,
        ], field.coordinate_system)
        return gradient
    raise ValueError(f"Unsupported coordinate system: {field.coordinate_system}")


def divergence_operator(field: VectorField) -> ScalarValue:
    field_space = field.apply_to_basis()
    # extend missing components with zeroes
    field_components = list(field_space.components) + [0] * (3 - len(field_space.components))
    if field.coordinate_system.coord_system_type == CoordinateSystem.System.CARTESIAN:
        x = field_space.coordinate_system.coord_system.base_scalars()[0]
        y = field_space.coordinate_system.coord_system.base_scalars()[1]
        z = field_space.coordinate_system.coord_system.base_scalars()[2]
        field_x = field_components[0]
        field_y = field_components[1]
        field_z = field_components[2]
        return diff(field_x, x) + diff(field_y, y) + diff(field_z, z)
    if field.coordinate_system.coord_system_type == CoordinateSystem.System.CYLINDRICAL:
        r = field_space.coordinate_system.coord_system.base_scalars()[0]
        theta = field_space.coordinate_system.coord_system.base_scalars()[1]
        z = field_space.coordinate_system.coord_system.base_scalars()[2]
        field_r = field_components[0]
        field_theta = field_components[1]
        field_z = field_components[2]
        return diff(field_r, r) + field_r / r + diff(field_theta, theta) / r + diff(field_z, z)
    if field.coordinate_system.coord_system_type == CoordinateSystem.System.SPHERICAL:
        r = field_space.coordinate_system.coord_system.base_scalars()[0]
        theta = field_space.coordinate_system.coord_system.base_scalars()[1]
        phi = field_space.coordinate_system.coord_system.base_scalars()[2]
        field_r = field_components[0]
        field_theta = field_components[1]
        field_phi = field_components[2]
        return diff(field_r, r) + 2 * field_r / r + diff(field_theta,
            theta) / (r * sin(phi)) + diff(field_phi, phi) / r + field_phi / (r * tan(phi))
    raise ValueError(f"Unsupported coordinate system: {field.coordinate_system}")


# Calculate Curl of the field, which is Cross(Nabla, Field)
def curl_operator(field: VectorField) -> VectorField:
    field_space = field.apply_to_basis()
    dimensions = 3
    if len(field_space.components) > dimensions:
        raise ValueError(
            f"Curl is only defined for {dimensions} dimensions. Got: {len(field_space.components)}")
    # extend missing components with zeroes
    field_components = list(field_space.components) + [0] * (3 - len(field_space.components))
    if field.coordinate_system.coord_system_type == CoordinateSystem.System.CARTESIAN:
        x = field_space.coordinate_system.coord_system.base_scalars()[0]
        y = field_space.coordinate_system.coord_system.base_scalars()[1]
        z = field_space.coordinate_system.coord_system.base_scalars()[2]
        field_x = field_components[0]
        field_y = field_components[1]
        field_z = field_components[2]
        field_rotor_vector = Vector([
            diff(field_z, y) - diff(field_y, z),
            diff(field_x, z) - diff(field_z, x),
            diff(field_y, x) - diff(field_x, y)
        ], field_space.coordinate_system)
        return VectorField.from_vector(field_rotor_vector)
    if field.coordinate_system.coord_system_type == CoordinateSystem.System.CYLINDRICAL:
        r = field_space.coordinate_system.coord_system.base_scalars()[0]
        theta = field_space.coordinate_system.coord_system.base_scalars()[1]
        z = field_space.coordinate_system.coord_system.base_scalars()[2]
        field_r = field_components[0]
        field_theta = field_components[1]
        field_z = field_components[2]
        field_rotor_vector = Vector([
            diff(field_z, theta) / r - diff(field_theta, z),
            diff(field_r, z) - diff(field_z, r),
            (diff(r * field_theta, r) - diff(field_r, theta)) / r
        ], field_space.coordinate_system)
        return VectorField.from_vector(field_rotor_vector)
    if field.coordinate_system.coord_system_type == CoordinateSystem.System.SPHERICAL:
        r = field_space.coordinate_system.coord_system.base_scalars()[0]
        theta = field_space.coordinate_system.coord_system.base_scalars()[1]
        phi = field_space.coordinate_system.coord_system.base_scalars()[2]
        field_r = field_components[0]
        field_theta = field_components[1]
        field_phi = field_components[2]
        field_rotor_vector = Vector([(diff(sin(phi) * field_theta, phi) - diff(field_phi, theta)) /
            (r * sin(phi)), (diff(r * field_phi, r) - diff(field_r, phi)) / r,
            (diff(field_r, theta) / sin(phi) - diff(r * field_theta, r)) / r],
            field_space.coordinate_system)
        return VectorField.from_vector(field_rotor_vector)
    raise ValueError(f"Unsupported coordinate system: {field.coordinate_system}")

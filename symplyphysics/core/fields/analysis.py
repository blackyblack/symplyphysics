from typing import Sequence
from sympy import Expr, integrate, simplify
from symplyphysics import Vector, dot_vectors, vector_magnitude, vector_unit
from ..dimensions import ScalarValue
from ..fields.operators import curl_operator, divergence_operator
from ..fields.vector_field import VectorField
from ..geometry.elements import parametrized_curve_element, parametrized_curve_element_magnitude, volume_element_magnitude
from ..geometry.normals import parametrized_curve_normal, parametrized_surface_normal
from ..fields.parameters import ParameterLimits


# trajectory should be array with projections to coordinates, eg [3 * cos(parameter), 3 * sin(parameter)]
def circulation_along_curve(field: VectorField, trajectory: Sequence[Expr],
    parameter_limits: ParameterLimits) -> ScalarValue:
    (parameter, parameter_from, parameter_to) = parameter_limits
    field_applied = field.apply(trajectory)
    curve_element_vector = parametrized_curve_element(Vector(trajectory, field.coordinate_system),
        parameter)
    integrand = dot_vectors(field_applied, curve_element_vector)
    circulation_value = integrate(integrand, (parameter, parameter_from, parameter_to))
    return simplify(circulation_value)


# calculate circulation along curve using surface that has this curve as a boundary
# surface should be array with projections to coordinates, eg [parameter1 * cos(parameter2), parameter1 * sin(parameter2)]
def circulation_along_surface_boundary(field: VectorField, surface: Sequence[Expr],
    parameter_and_limits1: ParameterLimits, parameter_and_limits2: ParameterLimits) -> ScalarValue:
    # circulation over surface is flux of curl of the field
    field_rotor_vector_field = curl_operator(field)
    return flux_across_surface(field_rotor_vector_field, surface, parameter_and_limits1,
        parameter_and_limits2)


# trajectory should be array with projections to coordinates, eg [3 * cos(parameter), 3 * sin(parameter)]
# trajectory and field should be 2-dimensional, on XY plane
def flux_across_curve(field: VectorField, trajectory: Sequence[Expr],
    parameter_limits: ParameterLimits) -> ScalarValue:
    if len(trajectory) > 2:
        raise ValueError(f"Trajectory should have at most 2 components, got {len(trajectory)}")
    (parameter, parameter_from, parameter_to) = parameter_limits
    field_applied = field.apply(trajectory)
    trajectory_vector = Vector(trajectory, field.coordinate_system)
    norm_vector = parametrized_curve_normal(trajectory_vector, parameter)
    norm_unit_vector = vector_unit(norm_vector)
    field_dot_norm_value = dot_vectors(field_applied, norm_unit_vector)
    curve_element_magnitude_value = parametrized_curve_element_magnitude(
        trajectory_vector, parameter)
    flux_value = integrate(field_dot_norm_value * curve_element_magnitude_value,
        (parameter, parameter_from, parameter_to))
    return simplify(flux_value)


# trajectory should be array with projections to coordinates, eg [3 * cos(parameter), 3 * sin(parameter)]
def flux_across_surface(field: VectorField, surface: Sequence[Expr],
    parameter_and_limits1: ParameterLimits, parameter_and_limits2: ParameterLimits) -> ScalarValue:
    (parameter1, parameter1_from, parameter1_to) = parameter_and_limits1
    (parameter2, parameter2_from, parameter2_to) = parameter_and_limits2
    # calculate SurfaceIntegral integrand, which is Dot(Field, dS)
    field_applied = field.apply(surface)
    surface_vector = Vector(surface, field.coordinate_system)
    surface_element_vector = parametrized_surface_normal(surface_vector, parameter1, parameter2)
    integrand = dot_vectors(field_applied, surface_element_vector)
    flux_value = integrate(integrand, (parameter1, parameter1_from, parameter1_to),
        (parameter2, parameter2_from, parameter2_to))
    return simplify(flux_value)


# flux across some curve, that is a surface boundary is double integral of divergence of the field
def flux_across_surface_boundary(field: VectorField, surface: Sequence[Expr],
    parameter_and_limits1: ParameterLimits, parameter_and_limits2: ParameterLimits) -> ScalarValue:
    (parameter1, parameter1_from, parameter1_to) = parameter_and_limits1
    (parameter2, parameter2_from, parameter2_to) = parameter_and_limits2
    field_divergence = divergence_operator(field)
    surface_vector = Vector(surface, field.coordinate_system)
    surface_element_vector = parametrized_surface_normal(surface_vector, parameter1, parameter2)
    surface_element_magnitude = vector_magnitude(surface_element_vector)
    flux_value = integrate(field_divergence * surface_element_magnitude,
        (parameter1, parameter1_from, parameter1_to), (parameter2, parameter2_from, parameter2_to))
    return simplify(flux_value)


# flux across some surface, that is a volume boundary is triple integral of divergence of the field
# over volume.
# Parametrized volumes are not supported. We define volume by the integral limits.
# Integration starts from the last limit, ie z_limits
def flux_across_volume_boundary(field: VectorField, x_limits: tuple[ScalarValue, ScalarValue],
    y_limits: tuple[ScalarValue, ScalarValue], z_limits: tuple[ScalarValue,
    ScalarValue]) -> ScalarValue:
    (x_from, x_to) = x_limits
    (y_from, y_to) = y_limits
    (z_from, z_to) = z_limits
    field_divergence = divergence_operator(field)
    x = field.coordinate_system.coord_system.base_scalars()[0]
    y = field.coordinate_system.coord_system.base_scalars()[1]
    z = field.coordinate_system.coord_system.base_scalars()[2]
    volume_element_magnitude_value = volume_element_magnitude(field.coordinate_system)
    flux_value = integrate(field_divergence * volume_element_magnitude_value, (z, z_from, z_to),
        (y, y_from, y_to), (x, x_from, x_to))
    return simplify(flux_value)

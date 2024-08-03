from sympy import diff, sympify
from symplyphysics import (QuantityVector, Symbol, add_cartesian_vectors, scale_vector, units,
    Quantity, validate_input, validate_output)
from symplyphysics.core.fields.operators import curl_operator
from symplyphysics.core.fields.vector_field import VectorField
from symplyphysics.core.dimensions import assert_equivalent_dimension
from symplyphysics.core.vectors.vectors import Vector

## Description
## The magnetic field circulation theorem states that an electric current and a change in electric induction generate
## a vortex magnetic field. Also known as Ampere's circuital law.

## Law is: curl(H(r, t)) = j + d(D(r, t)) / dt, where
## H - magnetic intensity (vector field),
## D - electric induction (vector field),
## j - conductivity current density vector,
## r - radius vector (i.e. x, y, z for point in cartesian coodrinates),
## t - time,
## d / dt - partial time derivative,
## curl - curl of vector field.

# Note:
## Curl operator is only defined in 3D space. Therefore magnetic intensity should not be defined
## as 4 dimensional field. Rather time should be a part of vector field component expression.

time = Symbol("time", units.time)


def conductivity_current_density_vector_law(magnetic_intensity: VectorField,
    electric_induction: VectorField) -> Vector:
    rotor_magnetic_intensity = curl_operator(magnetic_intensity)
    rotor_magnetic_intensity_space = rotor_magnetic_intensity.apply_to_basis()
    electric_induction_space = electric_induction.apply_to_basis()
    electric_induction_time_derivative = [
        diff(c, time) for c in electric_induction_space.components
    ]
    electric_induction_time_derivative_vector = Vector(electric_induction_time_derivative,
        electric_induction.coordinate_system)
    return add_cartesian_vectors(rotor_magnetic_intensity_space,
        scale_vector(-1, electric_induction_time_derivative_vector))


def magnetic_intensity_rotor_law(electric_induction: VectorField,
    conductivity_current_density: Vector) -> VectorField:
    electric_induction_space = electric_induction.apply_to_basis()
    electric_induction_time_derivative = [
        diff(c, time) for c in electric_induction_space.components
    ]
    electric_induction_time_derivative_vector = Vector(electric_induction_time_derivative,
        electric_induction.coordinate_system)
    sum_vectors = add_cartesian_vectors(conductivity_current_density,
        electric_induction_time_derivative_vector)
    return VectorField.from_vector(sum_vectors)


def electric_induction_time_derivative_law(magnetic_intensity: VectorField,
    conductivity_current_density: Vector) -> VectorField:
    rotor_magnetic_intensity = curl_operator(magnetic_intensity)
    rotor_magnetic_intensity_space = rotor_magnetic_intensity.apply_to_basis()
    sub_vectors = add_cartesian_vectors(rotor_magnetic_intensity_space,
        scale_vector(-1, conductivity_current_density))
    return VectorField.from_vector(sub_vectors)


@validate_input(cartesian_point_=units.length, time_=units.time)
@validate_output(units.current / units.area)
def calculate_conductivity_current_density_at_point(magnetic_intensity_: VectorField,
    electric_induction_: VectorField, cartesian_point_: tuple[Quantity, Quantity,
    Quantity], time_: Quantity) -> QuantityVector:
    magnetic_intensity_vector = magnetic_intensity_.apply(cartesian_point_)
    for i, c in enumerate(magnetic_intensity_vector.components):
        assert_equivalent_dimension(c, f"magnetic_intensity_vector[{i}]",
            "calculate_conductivity_current_density_at_point", units.current / units.length)

    electric_induction_vector = electric_induction_.apply(cartesian_point_)
    for i, c in enumerate(electric_induction_vector.components):
        expr = sympify(c)
        expr = expr.subs(time, time_)
        assert_equivalent_dimension(expr, f"electric_induction_vector[{i}]",
            "calculate_conductivity_current_density_at_point", units.charge / units.area)

    result_vector = conductivity_current_density_vector_law(magnetic_intensity_,
        electric_induction_)
    base_scalars = result_vector.coordinate_system.coord_system.base_scalars()
    return QuantityVector.from_base_vector(result_vector,
        subs={
        base_scalars[0]: cartesian_point_[0],
        base_scalars[1]: cartesian_point_[1],
        base_scalars[2]: cartesian_point_[2],
        time: time_
        })

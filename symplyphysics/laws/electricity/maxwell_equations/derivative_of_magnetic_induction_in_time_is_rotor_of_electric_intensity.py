from symplyphysics import (units, Symbol, Quantity, QuantityVector, validate_input, validate_output)
from symplyphysics.core.dimensions import assert_equivalent_dimension
from symplyphysics.core.fields.operators import curl_operator
from symplyphysics.core.fields.vector_field import VectorField
from symplyphysics.core.vectors.arithmetics import diff_cartesian_vector, scale_vector

## Description
## Faraday's law of induction states that a change in magnetic induction generates a vortex electric field.
## This law is valid for any magnetic field that changes over time.

## Law is: curl(E(r, t)) = -d(B(r, t)) / dt, where
## E - the field of electric intensity,
## B - the field of magnetic induction,
## r - the vector of a point in space,
## t - time,
## curl - curl of the vector of magnetic induction,
## d / dt - time partial derivative.

# Links
## Wikipedia, fourth line in table <https://en.wikipedia.org/wiki/Maxwell%27s_equations#Macroscopic_formulation>
## Physics LibreTexts, formula 15.7.1 <https://phys.libretexts.org/Bookshelves/Electricity_and_Magnetism/Electricity_and_Magnetism_(Tatum)/15%3A_Maxwell's_Equations/15.07%3A_Maxwell's_Fourth_Equation>

time = Symbol('time', units.time)


def magnetic_induction_derivative_law(electric_intensity_: VectorField) -> VectorField:
    return curl_operator(electric_intensity_)


def electric_intensity_curl_law(magnetic_induction_: VectorField) -> VectorField:
    magnetic_induction_space = magnetic_induction_.apply_to_basis()
    magnetic_induction_time_derivative = diff_cartesian_vector(magnetic_induction_space, time)
    magnetic_induction_time_derivative = scale_vector(-1, magnetic_induction_time_derivative)
    return VectorField.from_vector(magnetic_induction_time_derivative)


@validate_input(cartesian_point_=units.length, time_=units.time)
@validate_output(units.magnetic_density / units.time)
def calculate_magnetic_induction_derivative_at_point(electric_intensity_: VectorField,
    cartesian_point_: tuple[Quantity, Quantity, Quantity], time_: Quantity) -> QuantityVector:
    electric_intensity_vector = electric_intensity_.apply(cartesian_point_)
    for i, c in enumerate(electric_intensity_vector.components):
        assert_equivalent_dimension(c, f"electric_intensity_vector[{i}]",
            "calculate_magnetic_induction_derivative_at_point", units.force / units.charge)

    result_field = magnetic_induction_derivative_law(electric_intensity_)
    result_space = result_field.apply_to_basis()
    base_scalars = result_field.coordinate_system.coord_system.base_scalars()
    return QuantityVector.from_base_vector(result_space,
        subs={
        base_scalars[0]: cartesian_point_[0],
        base_scalars[1]: cartesian_point_[1],
        base_scalars[2]: cartesian_point_[2],
        time: time_
        })

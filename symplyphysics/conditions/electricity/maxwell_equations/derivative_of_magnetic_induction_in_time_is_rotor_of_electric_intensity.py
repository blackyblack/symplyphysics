from sympy import Eq, Derivative, sin, sympify, Matrix
from sympy.core import Equality
from symplyphysics import (units, Quantity, Symbol, validate_input)
from symplyphysics.core.fields.operators import curl_operator
from symplyphysics.core.fields.vector_field import VectorField
from symplyphysics.core.dimensions import assert_equivalent_dimension
from symplyphysics.core.vectors.vectors import Vector

## Description
## Faraday's law of induction states that a change in magnetic induction generates a vortex electric field.

## Law is: curl(E(r, t)) = -d(B(r, t)) / dt, where
## E - the field of electric intensity,
## B - the field of magnetic induction,
## r - the vector of a point in space,
## t - time,
## curl - curl of the vector of magnetic induction,
## d / dt - time partial derivative.

time = Symbol('time', units.time)


# def faradays_law_of_induction(electric_intensity_: VectorField,
#                                             magnetic_induction_: VectorField,
#                                             time_: Quantity) -> tuple[Equality, Equality, Equality]:
#     curl_operator_electric_intensity = curl_operator(electric_intensity_)

#     x = curl_operator_electric_intensity.coordinate_system.coord_system.base_scalars()[0]
#     y = curl_operator_electric_intensity.coordinate_system.coord_system.base_scalars()[1]
#     z = curl_operator_electric_intensity.coordinate_system.coord_system.base_scalars()[2]
#     dict_xyz = {
#         x: x_coor,
#         y: y_coor,
#         z: z_coor,
#         time_: time,
#     }

#     law_x = Eq(curl_operator_electric_intensity.apply_to_basis().components[0].subs(dict_xyz), -Derivative(magnetic_induction_.apply_to_basis().components[0].subs(dict_xyz), (time, 1)))
#     law_y = Eq(curl_operator_electric_intensity.apply_to_basis().components[1].subs(dict_xyz), -Derivative(magnetic_induction_.apply_to_basis().components[1].subs(dict_xyz), (time, 1)))
#     law_z = Eq(curl_operator_electric_intensity.apply_to_basis().components[2].subs(dict_xyz), -Derivative(magnetic_induction_.apply_to_basis().components[2].subs(dict_xyz), (time, 1)))
#     return law_x, law_y, law_z


def faradays_law_of_induction(electric_intensity_: VectorField, 
                              magnetic_induction_: VectorField) -> tuple[Equality, Equality, Equality]:
    curl_operator_electric_intensity = curl_operator(electric_intensity_)
    # curl_operator_electri_intensity_space = curl_operator_electric_intensity.apply_to_basis()

    magnetic_induction_space = magnetic_induction_.apply_to_basis()
    magnetic_induction_time_derivative = [-Derivative(c, time) for c in magnetic_induction_space.components]
    magnetic_induction_time_derivative_vector = Vector(magnetic_induction_time_derivative, magnetic_induction_.coordinate_system)
    
    # law_x = Eq(curl_operator_electri_intensity_space.components[0], -Derivative(magnetic_induction_space.components[0], time))
    # law_y = Eq(curl_operator_electri_intensity_space.components[1], -Derivative(magnetic_induction_space.components[1], time))
    # law_z = Eq(curl_operator_electri_intensity_space.components[2], -Derivative(magnetic_induction_space.components[2], time))
    return Eq(Matrix(curl_operator_electric_intensity.apply_to_basis().components), Matrix(magnetic_induction_time_derivative))
    return Eq(rotor_electric_intensity, VectorField.from_vector(magnetic_induction_time_derivative_vector))

@validate_input(cartesian_point_=units.length,
                time_=units.time)
def get_faradays_law_of_induction(electric_intensity_: VectorField,
                                                    magnetic_induction_: VectorField,
                                                    time_: Quantity,
                                                    cartesian_point_: tuple[Quantity, Quantity, Quantity]) -> tuple[Equality, Equality, Equality]:
    electric_intensity_vector = electric_intensity_.apply(cartesian_point_)
    for i, c in enumerate(electric_intensity_vector.components):
        expr = sympify(c)
        expr = expr.subs(time, time_)
        assert_equivalent_dimension(expr, f"electric_intensity_vector[{i}]",
            "get_faradays_law_of_induction", units.voltage / units.length)

    magnetic_induction_vector = magnetic_induction_.apply(cartesian_point_)
    for i, c in enumerate(magnetic_induction_vector.components):
        expr = sympify(c)
        expr = expr.subs(time, time_)
        assert_equivalent_dimension(expr, f"magnetic_induction_vector[{i}]",
            "get_faradays_law_of_induction", units.magnetic_density)

    result = faradays_law_of_induction(electric_intensity_, magnetic_induction_)

    return result


electric_amplitude = Quantity(10e3 * units.volt / units.meter)
mag_amplitude = Quantity(1e3 * units.tesla)
time__ = Quantity(1 * units.second)
x__ = Quantity(0.493 * units.meter)
y__ = Quantity(0.5 * units.meter)
z__ = Quantity(0.01249 * units.meter)
electric_intensity__ = VectorField(lambda point: [sin(point.z / Quantity(1 * units.meter)) * electric_amplitude * sin(time / Quantity(1 * units.second)),
                                                sin(point.x / Quantity(1 * units.meter)) * electric_amplitude * sin(time / Quantity(1 * units.second)),
                                                sin(point.y / Quantity(1 * units.meter)) * electric_amplitude * sin(time / Quantity(1 * units.second))])
magnetic_induction__ = VectorField(lambda point: [sin(point.x / Quantity(1 * units.meter)) * mag_amplitude * sin(time / Quantity(1 * units.second)),
                                                    sin(point.y / Quantity(1 * units.meter)) * mag_amplitude * sin(time / Quantity(1 * units.second)),
                                                    sin(point.z / Quantity(1 * units.meter)) * mag_amplitude * sin(time / Quantity(1 * units.second))])

laws = get_faradays_law_of_induction(electric_intensity__, magnetic_induction__, time__, (x__, y__, z__))
res = 1
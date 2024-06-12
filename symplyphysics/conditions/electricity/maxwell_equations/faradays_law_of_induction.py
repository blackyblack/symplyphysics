from sympy import Eq, Derivative
from sympy.core import Equality
from symplyphysics import (units, Quantity, Symbol, validate_input)
from symplyphysics.core.fields.operators import curl_operator as rotor
from symplyphysics.core.fields.vector_field import VectorField
from symplyphysics.core.dimensions import assert_equivalent_dimension

## Description
## Faraday's law of induction states that a change in magnetic induction generates a vortex electric field.

## Law is: rot(E(x, y, z, t)) = -d(B(x, y, z, t)) / dt, where
## E - the vector of electric intensity,
## B - the vector of magnetic induction,
## x, y, z - coordinates in a rectangular coordinate system,
## t - tiime,
## rotor - rotor of the vector of magnetic induction,
## d / dt - time derivative.

x_coor = Symbol('x_coor', units.length)
y_coor = Symbol('y_coor', units.length)
z_coor = Symbol('z_coor', units.length)
time = Symbol('time', units.time)


def faradays_law_of_induction(electric_intensity_: VectorField,
                                            magnetic_induction_: VectorField,
                                            time_: Quantity) -> tuple[Equality, Equality, Equality]:
    rotor_electric_intensity = rotor(electric_intensity_)

    x = rotor_electric_intensity.coordinate_system.coord_system.base_scalars()[0]
    y = rotor_electric_intensity.coordinate_system.coord_system.base_scalars()[1]
    z = rotor_electric_intensity.coordinate_system.coord_system.base_scalars()[2]
    dict_xyz = {
        x: x_coor,
        y: y_coor,
        z: z_coor,
        time_: time,
    }

    law_x = Eq(rotor_electric_intensity.apply_to_basis().components[0].subs(dict_xyz), -Derivative(magnetic_induction_.apply_to_basis().components[0].subs(dict_xyz), (time, 1)))
    law_y = Eq(rotor_electric_intensity.apply_to_basis().components[1].subs(dict_xyz), -Derivative(magnetic_induction_.apply_to_basis().components[1].subs(dict_xyz), (time, 1)))
    law_z = Eq(rotor_electric_intensity.apply_to_basis().components[2].subs(dict_xyz), -Derivative(magnetic_induction_.apply_to_basis().components[2].subs(dict_xyz), (time, 1)))
    return law_x, law_y, law_z


@validate_input(cartesian_point_=units.length,
                time_=units.time)
def get_faradays_law_of_induction(electric_intensity_: VectorField,
                                                    magnetic_induction_: VectorField,
                                                    time_: Quantity,
                                                    cartesian_point_: tuple[Quantity, Quantity, Quantity]) -> tuple[Equality, Equality, Equality]:
    electric_intensity_vector = electric_intensity_.apply(cartesian_point_)
    for i, c in enumerate(electric_intensity_vector.components):
        assert_equivalent_dimension(c, f"electric_intensity_vector[{i}]",
            "get_faradays_law_of_induction", units.voltage / units.length)

    magnetic_induction_vector = magnetic_induction_.apply(cartesian_point_)
    for i, c in enumerate(magnetic_induction_vector.components):
        assert_equivalent_dimension(c, f"magnetic_induction_vector[{i}]",
            "get_faradays_law_of_induction", units.magnetic_density)

    result = faradays_law_of_induction(electric_intensity_, magnetic_induction_, time_)

    return result

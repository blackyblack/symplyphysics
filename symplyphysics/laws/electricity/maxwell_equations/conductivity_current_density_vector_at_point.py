from sympy import diff
from symplyphysics import (units, Quantity, validate_input)
from symplyphysics.core.fields.operators import curl_operator as rotor
from symplyphysics.core.fields.vector_field import VectorField
from symplyphysics.core.dimensions import assert_equivalent_dimension
from symplyphysics.core.vectors.vectors import Vector

## Description
## The magnetic field circulation theorem states that an electric current and a change in electrical induction generate
## a vortex magnetic field.

## Law is: rot(H(x, y, z, t)) = j + d(D(x, y, z, t)) / dt, where
## H - the vector of magnetic intensity,
## D - the vector of electrical induction,
## j - conductivity current density vector,
## x, y, z - coordinates in a rectangular coordinate system,
## t - time,
## d / dt - time derivative,
## rot - rotor of vector.


def conductivity_current_density_vector_law(magnetic_intensity: VectorField,
                                            electric_induction: VectorField,
                                            cartesian_point_: tuple[Quantity, Quantity, Quantity],
                                            time: Quantity) -> Vector:
    rotor_magnetic_intensity = rotor(magnetic_intensity)
    rotor_magnetic_intensity = rotor_magnetic_intensity.apply(cartesian_point_)

    electric_induction_vector = electric_induction.apply(cartesian_point_)

    current_vector = Vector([Quantity(rotor_magnetic_intensity.components[0] - diff(electric_induction_vector.components[0], time)),
                             Quantity(rotor_magnetic_intensity.components[1] - diff(electric_induction_vector.components[1], time)),
                             Quantity(rotor_magnetic_intensity.components[2] - diff(electric_induction_vector.components[2], time))])
    return current_vector
    

@validate_input(cartesian_point_=units.length,
                time_=units.time)
def calculate_conductivity_current_density_at_point(magnetic_intensity_: VectorField,
                                                    electric_induction_: VectorField,
                                                    cartesian_point_: tuple[Quantity, Quantity, Quantity],
                                                    time_: Quantity) -> Vector:
    magnetic_intensity_vector = magnetic_intensity_.apply(cartesian_point_)
    for i, c in enumerate(magnetic_intensity_vector.components):
        assert_equivalent_dimension(c, f"magnetic_intensity_vector[{i}]",
            "calculate_conductivity_current_density_at_point", units.current / units.length)


    electric_induction_vector = electric_induction_.apply(cartesian_point_)
    for i, c in enumerate(electric_induction_vector.components):
        assert_equivalent_dimension(c, f"electric_induction_vector[{i}]",
            "calculate_conductivity_current_density_at_point", units.charge / units.area)

    result = conductivity_current_density_vector_law(magnetic_intensity_, electric_induction_, cartesian_point_, time_)
    for i, c in enumerate(result.components):
        assert_equivalent_dimension(c, f"result[{i}]",
            "calculate_conductivity_current_density_at_point", units.current / units.area)

    return result



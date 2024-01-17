from symplyphysics import (
    Symbol,
    units,
    validate_input,
    validate_output,
    Quantity,
    scale_vector,
    list_of_quantities,
)
from symplyphysics.core.vectors.vectors import QuantityVector, Vector

# Description
## The dipole moment of a system of two point charges, one positive (q) and one negative (-q) is collinear with the
## displacement vector pointing from the negative charge to the positive charge. For a further explanation
## of what a dipole moment means, see (electric_dipole_moment.py)[../electric_dipole_moment.py].

# Law: p = q * l
## p - vector of dipole moment
## q - electric charge
## l - displacement vector from the negative charge to the positive charge

charge = Symbol("charge", units.charge)


def dipole_moment_law(displacement_vector_: Vector) -> Vector:
    return scale_vector(charge, displacement_vector_)


@validate_input(charge_=charge, displacement_vector_=units.length)
@validate_output(units.charge * units.length)
def calculate_dipole_moment(charge_: Quantity, displacement_vector_: QuantityVector) -> QuantityVector:
    result_dipole_moment = dipole_moment_law(displacement_vector_)
    dipole_moment_components = list_of_quantities(
        result_dipole_moment.components,
        {charge: charge_},
    )
    return QuantityVector(dipole_moment_components, displacement_vector_.coordinate_system)

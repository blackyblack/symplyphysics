"""
Electric dipole moment is charge times displacement
===================================================

The vector of electric dipole moment is a vector whose magnitude describes the
:doc:`electric dipole moment <laws.electricity.electric_dipole_moment_is_charge_times_distance>`
of the system. It is collinear to the vector connecting the two point charges.

**Conditions:**

#. The charges must be equal by magnitude and have opposite signs.
"""

from symplyphysics import (
    Symbol,
    units,
    validate_input,
    validate_output,
    Quantity,
    scale_vector,
    QuantityVector,
    Vector,
)

charge = Symbol("charge", units.charge)
"""
Magnitude of the electric charge of the point charges.

Symbol:
    :code:`q`
"""


def dipole_moment_law(displacement_vector_: Vector) -> Vector:
    r"""
    Vector of electric dipole moment via displacement vector.

    Law:
        :code:`p = q * d`

    Latex:
        .. math::
            \vec p = q \vec d

    :param displacement\_vector\_: vector pointing from the negative charge to the positive charge

        Symbol: :code:`d`

        Latex: :math:`\vec d`

        Dimension: *length*

    :return: vector of electric dipole moment

        Symbol: :code:`p`

        Latex: :math:`\vec p`

        Dimension: *charge* * *length*
    """

    return scale_vector(charge, displacement_vector_)


@validate_input(charge_=charge, displacement_vector_=units.length)
@validate_output(units.charge * units.length)
def calculate_dipole_moment(charge_: Quantity,
    displacement_vector_: QuantityVector) -> QuantityVector:
    result_vector = dipole_moment_law(displacement_vector_.to_base_vector())
    return QuantityVector.from_base_vector(result_vector, subs={charge: charge_})

"""
Electric field is force over test charge (Vector)
=================================================

Electric field at a point in space can be found by placing there a test charge and measuring
the electrostatic force that is applied to it.
"""

from symplyphysics import (
    Quantity,
    QuantityVector,
    scale_vector,
    Symbol,
    units,
    validate_input,
    validate_output,
    Vector,
)

test_charge = Symbol("test_charge", units.charge)
r"""
Value of the test charge.

Symbol:
    :code:`q0`

Latex:
    :math:`q_0`
"""


def electric_field_law(electrostatic_force_: Vector) -> Vector:
    r"""
    Electric field via electrostatic force.

    Law:
        :code:`E = F / q0`

    Latex:
        .. math::
            \vec E = \frac{\vec F}{q_0}

    :param electrostatic\_force\_: vector of electrostatic force

        Symbol: :code:`F`

        Latex: :math:`\vec F`

        Dimension: *force*

    :return: vector of electric field

        Symbol: :code:`E`

        Latex: :math:`\vec E`

        Dimension: *force* / *charge*
    """

    return scale_vector(1 / test_charge, electrostatic_force_)


def electrostatic_force_law(electric_field_: Vector) -> Vector:
    r"""
    Electrostatic force via electric field.

    Law:
        :code:`F = q0 * E`

    Latex:
        .. math::
            \vec F = q_0 \vec E

    :param electric\_field\_: vector of electric field

        Symbol: :code:`E`

        Latex: :math:`\vec E`

        Dimension: *force* / *charge*

    :return: vector of electrostatic force

        Symbol: :code:`F`

        Latex: :math:`\vec F`

        Dimension: *force*
    """

    return scale_vector(test_charge, electric_field_)


@validate_input(electrostatic_force_=units.force, test_charge_=test_charge)
@validate_output(units.force / units.charge)
def calculate_electric_field(electrostatic_force_: QuantityVector,
    test_charge_: Quantity) -> QuantityVector:
    result_vector = electric_field_law(electrostatic_force_.to_base_vector())
    return QuantityVector.from_base_vector(result_vector, subs={test_charge: test_charge_})


@validate_input(electric_field_=units.force / units.charge, test_charge_=test_charge)
@validate_output(units.force)
def calculate_electrostatic_force(electric_field_: QuantityVector,
    test_charge_: Quantity) -> QuantityVector:
    result_vector = electrostatic_force_law(electric_field_.to_base_vector())
    return QuantityVector.from_base_vector(result_vector, subs={test_charge: test_charge_})

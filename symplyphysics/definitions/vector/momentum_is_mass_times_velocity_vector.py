"""
Momentum is mass times velocity (Vector)
========================================

An object's *linear momentum* is a vector quantity defined as the product of its mass and velocity vector.
"""

from symplyphysics import (
    units,
    symbols,
    Vector,
    Quantity,
    QuantityVector,
    scale_vector,
    validate_input,
    validate_output,
)

mass = symbols.basic.mass
"""
:attr:`~symplyphysics.symbols.basic.mass` of the object.
"""


def momentum_definition(velocity_: Vector) -> Vector:
    r"""
    Vector of *linear momentum*.

    Law:
        :code:`p = m * v`

    Latex:
        .. math::
            \vec p = m \vec v

    :param velocity\_: vector of velocity of the object.

        Symbol: :code:`v`

        Latex: :math:`\vec v`

        Dimension: *velocity*

    :return: vector of linear momentum.

        Symbol: :code:`p`

        Latex: :math:`\vec p`

        Dimension: *momentum*
    """

    return scale_vector(mass, velocity_)


def velocity_law(momentum_: Vector) -> Vector:
    r"""
    Vector of *velocity*

    Law:
        :code:`v = p / m`

    Latex:
        .. math::
            \vec v = \frac{\vec p}{m}

    :param momentum\_: vector of linear momentum.

        Symbol: :code:`p`

        Latex: :math:`\vec p`

        Dimension: *momentum*

    :return: vector of velocity of the object.

        Symbol: :code:`v`

        Latex: :math:`\vec v`

        Dimension: *velocity*
    """

    return scale_vector(1 / mass, momentum_)


@validate_input(mass_=mass, velocity_=units.velocity)
@validate_output(units.momentum)
def calculate_momentum(mass_: Quantity, velocity_: QuantityVector) -> QuantityVector:
    result_vector = momentum_definition(velocity_.to_base_vector())
    return QuantityVector.from_base_vector(result_vector, subs={mass: mass_})


@validate_input(mass_=mass, momentum_=units.momentum)
@validate_output(units.velocity)
def calculate_velocity(mass_: Quantity, momentum_: QuantityVector) -> QuantityVector:
    result_vector = velocity_law(momentum_.to_base_vector())
    return QuantityVector.from_base_vector(result_vector, subs={mass: mass_})

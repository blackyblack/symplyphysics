r"""
Kinetic energy via angular momentum and angular velocity
========================================================

Kinetic energy of a body rotating around a fixed or instantaneous axis depends on its
angular momentum and angular velocity.
"""

from sympy import Expr
from symplyphysics import (
    units,
    angle_type,
    Quantity,
    QuantityVector,
    Vector,
    dot_vectors,
    validate_input,
    validate_output,
)


def kinetic_energy_law(
    angular_momentum_: Vector,
    angular_velocity_: Vector,
) -> Expr:
    r"""
    Kinetic energy of a rotating body.

    **Notation:**

    #. :math:`\left(\vec a, \vec b)` (:code:`dot(a, b)`) is the dot product between vectors
       :math:`\vec a` and :math:`\vec b`.

    Law:
        :code:`K = 1/2 * dot(L, w)`

    Latex:
        :math:`K = \frac{1}{2} \left(\vec L, \vec \omega \right)`

    :param angular_momentum\_: angular momentum of the body.

        Symbol: :code:`L`

        Latex: :math:`\vec L`

        Dimension: *length* * *momentum*

    :param angular_velocity\_: angular velocity of the body.

        Symbol: :code:`w`

        Latex: :math:`\vec \omega`

        Dimension: *angle* / *time*

    :return: kinetic energy of the body

        Symbol: :code:`K`

        Dimension: *energy*
    """

    return dot_vectors(angular_momentum_, angular_velocity_) / 2


@validate_input(
    angular_momentum_=units.length * units.momentum,
    angular_velocity_=angle_type / units.time,
)
@validate_output(units.energy)
def calculate_kinetic_energy(
    angular_momentum_: QuantityVector,
    angular_velocity_: QuantityVector,
) -> Quantity:
    result = kinetic_energy_law(
        angular_momentum_.to_base_vector(),
        angular_velocity_.to_base_vector(),
    )
    return Quantity(result)

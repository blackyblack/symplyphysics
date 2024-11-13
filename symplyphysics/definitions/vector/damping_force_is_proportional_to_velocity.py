"""
Damping force is proportional to velocity
=========================================

Damping force is an external (relative to an object) force that drains energy from the object,
reducing the motion of the object. It is a model used, for example, to describe the motion
of an oscillator.
"""

from symplyphysics import (
    units,
    Symbol,
    Quantity,
    Vector,
    QuantityVector,
    scale_vector,
    validate_input,
    validate_output,
)

damping_constant = Symbol("damping_constant", units.mass / units.time)
"""
Non-negative damping constant.

Symbol:
    :code:`b`
"""


def damping_force_definition(velocity_: Vector) -> Vector:
    r"""
    Vector of *damping force* exerted on the object.

    Law:
        :code:`F = -1 * b * v`

    Latex:
        .. math::
            \vec F = -b \vec v

    :param velocity\_: velocity vector of the object.

        Symbol: :code:`v`

        Latex: :math:`\vec v`

        Dimension: *velocity*

    :return: vector of damping :symbols:`force`.

        Symbol: :code:`F`

        Latex: :math:`\vec F`

        Dimension: *force*
    """

    return scale_vector(-1 * damping_constant, velocity_)


def velocity_law(damping_force_: Vector) -> Vector:
    r"""
    *Velocity* of the object which the damping force is exerted on.
    
    Law:
        :code:`v = -1/b * F`

    Latex:
        .. math::
            \vec v = - \frac{\vec F}{b}
    
    :param damping_force\_: damping :symbols:`force` exerted on the object.

        Symbol: :code:`F`

        Latex: :math:`\vec F`

        Dimension: *force*

    :return: velocity vector of the object.

        Symbol: :code:`v`

        Latex: :math:`\vec v`

        Dimension: *velocity*
    """

    return scale_vector(-1 / damping_constant, damping_force_)


@validate_input(damping_constant_=damping_constant, velocity_=units.velocity)
@validate_output(units.force)
def calculate_damping_force(damping_constant_: Quantity,
    velocity_: QuantityVector) -> QuantityVector:
    result_vector = damping_force_definition(velocity_.to_base_vector())
    return QuantityVector.from_base_vector(result_vector,
        subs={damping_constant: damping_constant_})


@validate_input(damping_constant_=damping_constant, damping_force_=units.force)
@validate_output(units.velocity)
def calculate_velocity(damping_constant_: Quantity,
    damping_force_: QuantityVector) -> QuantityVector:
    result_vector = velocity_law(damping_force_.to_base_vector())
    return QuantityVector.from_base_vector(result_vector,
        subs={damping_constant: damping_constant_})

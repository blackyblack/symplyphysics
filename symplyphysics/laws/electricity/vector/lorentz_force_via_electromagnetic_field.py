r"""
Lorentz force via electromagnetic field
=======================================

The **Lorentz force law** states that a charged particle moving in an electromagnetic
field experiences a force that depends on the values of the electric field and the
magnetic field.

**Notation:**

#. :math:`\vec a \times \vec b` (:code:`cross(a, b)`) is the cross product between
   :math:`\vec a` and :math:`\vec b`.
#. :math:`\lVert \vec a \rVert` (:code:`norm(a)`) is the Euclidean norm of :math:`\vec a`.
#. :math:`|x|` (:code:`abs(x)`) is the absolute value of :math:`x`.

**Notes:**

#. This law is valid even in the relativistic case.
#. This law works only in principle because a real particle would generate its own
   electromagnetic field that would interact with the external one which would alter
   the electromagnetic force it experiences.
"""

from sympy import Expr
from symplyphysics import (
    Quantity,
    QuantityVector,
    scale_vector,
    add_cartesian_vectors,
    subtract_cartesian_vectors,
    cross_cartesian_vectors,
    vector_magnitude,
    Symbol,
    units,
    validate_input,
    validate_output,
    Vector,
)

charge = Symbol("charge", units.charge)
"""
Value of the electric charge of the test particle.

Symbol:
    :code:`q`
"""


def lorentz_force_law(
    electric_field_: Vector,
    magnetic_flux_density_: Vector,
    velocity_: Vector,
) -> Vector:
    r"""
    Lorentz force via electric and magnetic fields, and velocity.

    Law:
        :code:`F = q * (E + cross(v, B))`

    Latex:
        .. math::
            \vec F = q \left( \vec E + \vec v \times \vec B \right)

    :param electric\_field\_: vector of electric field

        Symbol: :code:`E`

        Latex: :math:`\vec E`

        Dimension: *voltage* / *length*
    
    :param magnetic\_flux\_density\_: pseudovector of magnetic flux density

        Symbol: :code:`B`

        Latex: :math:`\vec B`

        Dimension: *magnetic flux density*

    :param velocity\_: vector of particle's velocity

        Symbol: :code:`v`

        Latex: :math:`\vec v`

        Dimension: *velocity*
    
    :return: Lorentz force acting on the charged particle

        Symbol: :code:`F`

        Latex: :math:`\vec F`

        Dimension: *force*
    """

    velocity_cross_magnetic_field_ = cross_cartesian_vectors(
        velocity_,
        magnetic_flux_density_,
    )

    total_field_ = add_cartesian_vectors(
        electric_field_,
        velocity_cross_magnetic_field_,
    )

    return scale_vector(charge, total_field_)


def electric_field_law(
    lorentz_force_: Vector,
    magnetic_flux_density_: Vector,
    velocity_: Vector,
) -> Vector:
    r"""
    Electric field via Lorentz force, magnetic field, and velocity.

    Law:
        :code:`E = F / q - cross(v, B)`

    Latex:
        .. math::
            \vec E = \frac{\vec F}{q} - \vec v \times \vec B

    :param lorentz\_force\_: Lorentz force acting on the charged particle

        Symbol: :code:`F`

        Latex: :math:`\vec F`

        Dimension: *force*

    :param magnetic\_flux\_density\_: pseudovector of magnetic flux density

        Symbol: :code:`B`

        Latex: :math:`\vec B`

        Dimension: *magnetic flux density*

    :param velocity\_: vector of particle's velocity

        Symbol: :code:`v`

        Latex: :math:`\vec v`

        Dimension: *velocity*

    :return: vector of electric field

        Symbol: :code:`E`

        Latex: :math:`\vec E`

        Dimension: *voltage* / *length*
    """

    velocity_cross_magnetic_field_ = cross_cartesian_vectors(
        velocity_,
        magnetic_flux_density_,
    )

    total_field_ = scale_vector(1 / charge, lorentz_force_)

    return subtract_cartesian_vectors(
        total_field_,
        velocity_cross_magnetic_field_,
    )


def charge_law(
    lorentz_force_: Vector,
    electric_field_: Vector,
    magnetic_flux_density_: Vector,
    velocity_: Vector,
) -> Expr:
    r"""
    Magnitude of the particle's charge via force and electromagnetic field.

    Law:
        :code:`abs(q) = norm(F) / norm(E + cross(v, B))`

    Latex:
        .. math::
            |q| = \frac{\lVert \vec F \rVert}{\left \lVert \vec E + \vec v \times \vec B \right \rVert}

    :param lorentz\_force\_: Lorentz force acting on the charged particle

        Symbol: :code:`F`

        Latex: :math:`\vec F`

        Dimension: *force*

    :param electric\_field\_: vector of electric field

        Symbol: :code:`E`

        Latex: :math:`\vec E`

        Dimension: *voltage* / *length*
    
    :param magnetic\_flux\_density\_: pseudovector of magnetic flux density

        Symbol: :code:`B`

        Latex: :math:`\vec B`

        Dimension: *magnetic flux density*

    :param velocity\_: vector of particle's velocity

        Symbol: :code:`v`

        Latex: :math:`\vec v`

        Dimension: *velocity*

    :return: magnitude of the test charge

        Symbol: :code:`q`

        Dimension: *charge*
    """

    velocity_cross_magnetic_field_ = cross_cartesian_vectors(
        velocity_,
        magnetic_flux_density_,
    )

    total_field_ = add_cartesian_vectors(
        electric_field_,
        velocity_cross_magnetic_field_,
    )

    return vector_magnitude(lorentz_force_) / vector_magnitude(total_field_)


@validate_input(
    charge_=charge,
    electric_field_=units.voltage / units.length,
    magnetic_flux_density_=units.magnetic_flux_density,
    velocity_=units.velocity,
)
@validate_output(units.force)
def calculate_lorentz_force(
    charge_: Quantity,
    electric_field_: QuantityVector,
    magnetic_flux_density_: QuantityVector,
    velocity_: QuantityVector,
) -> QuantityVector:
    result = lorentz_force_law(
        electric_field_=electric_field_.to_base_vector(),
        magnetic_flux_density_=magnetic_flux_density_.to_base_vector(),
        velocity_=velocity_.to_base_vector(),
    )
    return QuantityVector.from_base_vector(
        result,
        subs={charge: charge_},
    )

"""
Acceleration from force and velocity
====================================

In special relativity, the Newton's second law does not hold in the classical form
:math:`\\vec F = m \\vec a`, but acceleration can still be expressed via force and velocity.

**Notation:**

#. :quantity_notation:`speed_of_light`.

**Conditions:**

#. This law applies to special relativity.

**Links:**

#. `Wikipedia, see paragraph <https://en.wikipedia.org/wiki/Acceleration_(special_relativity)#Acceleration_and_force>`__.
"""

from sympy import Eq

from symplyphysics import Quantity, validate_input, validate_output, symbols, quantities
from symplyphysics.definitions import lorentz_factor as lorentz_factor_def

from symplyphysics.core.experimental.vectors import clone_as_vector_symbol, VectorDot, VectorNorm
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

acceleration = clone_as_vector_symbol(symbols.acceleration)
"""
Vector of the body's :symbols:`acceleration`.
"""

rest_mass = symbols.rest_mass
"""
:symbols:`rest_mass` of the body.
"""

force = clone_as_vector_symbol(symbols.force)
"""
Vector of the :symbols:`force` exerted on the body.
"""

velocity = clone_as_vector_symbol(symbols.speed)
"""
Vector of the body's velocity. See :symbols:`speed`.
"""

lorentz_factor = symbols.lorentz_factor
"""
:symbols:`lorentz_factor`.
"""

law = Eq(
    acceleration,
    (rest_mass * lorentz_factor)**(-1) *
    (force - quantities.speed_of_light**(-2) * VectorDot(force, velocity) * velocity),
)
"""
:laws:symbol::

:laws:latex::
"""

# TODO: prove that this is derivable from `./force_acceleration_relation.py` and vice versa


@validate_input(rest_mass_=rest_mass, force_=force, velocity_=velocity)
@validate_output(acceleration)
def calculate_acceleration(
    rest_mass_: Quantity,
    force_: QuantityCoordinateVector,
    velocity_: QuantityCoordinateVector,
) -> QuantityCoordinateVector:
    lorentz_factor_ = lorentz_factor_def.definition.rhs.subs(
        lorentz_factor_def.speed,
        VectorNorm(velocity_),
    )

    result = law.rhs.subs({
        rest_mass: rest_mass_,
        lorentz_factor: lorentz_factor_,
        force: force_,
        velocity: velocity_,
    })

    return QuantityCoordinateVector.from_expr(result)

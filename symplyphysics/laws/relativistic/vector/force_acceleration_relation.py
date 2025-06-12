"""
Force from acceleration and velocity
====================================

In special relativity, the Newton's second law does not hold in the classical form
:math:`\\vec F = m \\vec a`, but force can still be expressed via acceleration and velocity.

**Conditions:**

#. This law applies to special relativity.

**Links:**

#. `Wikipedia, see paragraph <https://en.wikipedia.org/wiki/Acceleration_(special_relativity)#Acceleration_and_force>`__.

..
    TODO: rename file
"""

from sympy import Eq, evaluate, root
from symplyphysics import Quantity, validate_input, validate_output, symbols, quantities

from symplyphysics.core.experimental.vectors import clone_as_vector_symbol, VectorDot
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

force = clone_as_vector_symbol(symbols.force)
"""
Vector of the :symbols:`force` exerted on the body.
"""

rest_mass = symbols.rest_mass
"""
:symbols:`rest_mass` of the body.
"""

tangential_acceleration = clone_as_vector_symbol(
    symbols.acceleration,
    display_symbol="a_t",
    display_latex="{\\vec a}_\\text{t}",
)
"""
Vector of the body's :symbols:`acceleration` tangential to the :attr:`~velocity` vector.
"""

normal_acceleration = clone_as_vector_symbol(
    symbols.acceleration,
    display_symbol="a_n",
    display_latex="{\\vec a}_\\text{n}",
)
"""
Vector of the body's :symbols:`acceleration` normal to the :attr:`~velocity` vector.
"""

velocity = clone_as_vector_symbol(symbols.speed)
"""
Vector of the body's velocity. See :symbols:`speed`.
"""

lorentz_factor = symbols.lorentz_factor
"""
:symbols:`lorentz_factor`.
"""

with evaluate(False):
    _tangential_force = lorentz_factor**3 * rest_mass * tangential_acceleration
    _normal_force = lorentz_factor * rest_mass * normal_acceleration

law = Eq(force, _tangential_force + _normal_force)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    rest_mass_=rest_mass,
    acceleration_=symbols.acceleration,
    velocity_=velocity,
)
@validate_output(force)
def calculate_force(
    rest_mass_: Quantity,
    acceleration_: QuantityCoordinateVector,
    velocity_: QuantityCoordinateVector,
) -> QuantityCoordinateVector:
    v_dot_v = VectorDot(velocity_, velocity_)

    tangential_acceleration_ = VectorDot(acceleration_, velocity_) / v_dot_v * velocity_

    normal_acceleration_ = acceleration_ - tangential_acceleration_

    # refer to `symplyphysics.definitions.lorentz_factor`
    # NOTE: make a law for calculating the Lorentz factor using a velocity vector
    lorentz_factor_ = root(1 - v_dot_v / quantities.speed_of_light**2, -2)

    result = law.rhs.subs({
        lorentz_factor: lorentz_factor_,
        rest_mass: rest_mass_,
        tangential_acceleration: tangential_acceleration_,
        normal_acceleration: normal_acceleration_,
    })

    return QuantityCoordinateVector.from_expr(result)

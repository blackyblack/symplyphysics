"""
Relativistic mass moment
========================

Mass moment is an additive physical quantity useful for deriving the Lorentz transformation of
angular momentum. For isolated systems, it is conserved in time, but unlike angular momentum, the
vector of mass moment is a polar ("ordinary") vector and therefore invariant under inversion.

**Links:**

#. `Wikipedia, end of paragraph <https://en.wikipedia.org/wiki/Relativistic_angular_momentum#Dynamic_mass_moment>`__.
"""

from sympy import Eq

from symplyphysics import Quantity, validate_input, validate_output, symbols, units
from symplyphysics.definitions import lorentz_factor as lorentz_factor_def

from symplyphysics.core.experimental.vectors import VectorSymbol, clone_as_vector_symbol, VectorNorm
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

mass_moment = VectorSymbol("N", units.mass * units.length)
"""
Vector of mass moment.
"""

rest_mass = symbols.rest_mass
"""
:symbols:`rest_mass` of the body.
"""

time = symbols.time
"""
:symbols:`time`.
"""

lorentz_factor = symbols.lorentz_factor
"""
:symbols:`lorentz_factor`.
"""

position_vector = clone_as_vector_symbol(symbols.distance_to_origin)
"""
Vector of the body's position. See :symbols:`distance_to_origin`.
"""

velocity = clone_as_vector_symbol(symbols.speed)
"""
Vector of the body's velocity. See :symbols:`speed`.
"""

law = Eq(
    mass_moment,
    rest_mass * lorentz_factor**2 * (position_vector - velocity * time),
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    rest_mass_=rest_mass,
    position_=position_vector,
    velocity_=velocity,
    time_=time,
)
@validate_output(mass_moment)
def calculate_mass_moment(
    rest_mass_: Quantity,
    position_: QuantityCoordinateVector,
    velocity_: QuantityCoordinateVector,
    time_: Quantity,
) -> QuantityCoordinateVector:
    lorentz_factor_ = lorentz_factor_def.definition.rhs.subs(
        lorentz_factor_def.speed,
        VectorNorm(velocity_),
    )

    result = law.rhs.subs({
        rest_mass: rest_mass_,
        lorentz_factor: lorentz_factor_,
        position_vector: position_,
        velocity: velocity_,
        time: time_,
    })

    return QuantityCoordinateVector.from_expr(result)


# UNIQUE_LAW_ID: 709
